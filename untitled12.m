%% cylinderConvergence_meanArea_exportSTL_REG_LOCALREFINE_FAST.m
% -------------------------------------------------------------------------
% Dynamic mesh convergence for meshed CYLINDERS (LOCAL refinement only)
%
% Uses YOUR method (per refinement level):
%   - Run ALL frames on current mesh
%   - Accumulate UNION of contact faces across frames
%   - Grow mask by N rings
%   - Refine ONLY that union region (midpoint split)
%   - Track OVERALL MEAN contact area per level
%   - Converge when mean %diff < 5% for 3 consecutive levels
%   - Export converged FEMUR mesh as STL (LOCAL coordinates)
%
% Speed + stability upgrades vs your current script:
%   (1) PRECOMPUTE master face normals + KDTree + bbox (ROI contact speed)
%   (2) FAST mask growth using FACE adjacency (shared edges) rather than rebuilding v2f every time
%   (3) OPTIONAL animation (animateEveryLevel=false is massively faster)
%   (4) FIXED alignment pipeline: flexion correction + CONSTANT registration (frame 1)
%       so you are NOT relying on raw tracking axes lining up with STL axes.
%
% -------------------------------------------------------------------------
clear; clc; close all;

%% ---------------- PATHS ----------------
scriptDir    = fileparts(mfilename('fullpath'));
projectRoot  = scriptDir;
addpath(genpath(fullfile(projectRoot,'src')));

dataDir  = fullfile(projectRoot,'data');
modelDir = fullfile(projectRoot,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');

assert(exist(matPath,'file')==2,  'Missing cylinders.mat');
assert(exist(femurSTL,'file')==2, 'Missing femur STL');
assert(exist(tibiaSTL,'file')==2, 'Missing tibia STL');

%% ---------------- OPTIONS ----------------
trialIndex = 1;
tol        = 0.20;      % mm
maxLevels  = 25;
minLevels  = 3;
convTolPct = 5;
convConsec = 3;

expandRings = 1;        % grow union contact region by N face-adjacency rings
framePause  = 0.02;

% BIG speed lever:
animateEveryLevel = false;   % true = show animation for every refinement level (SLOW)
animateFinalOnly  = true;    % if animateEveryLevel=false, animate only the converged level
debugFrame1Align  = true;    % show frame 1 tibia+femur overlay after registration

% Contact function preference:
preferROI = true;           % uses computeContactArea_STS_hybrid_ROI if present

%% ---------------- LOAD DATA ----------------
S  = load(matPath);
tr = S.trajectories(trialIndex);
fTt = tr.Transform.fTt;      % tibia in femur
nFrames = size(fTt,3);

%% ---------------- LOAD STLs ----------------
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vmov0 = Sf.Points;  Fmov0 = Sf.ConnectivityList;   % femur (moving mesh in STL coords)
Vfix  = St.Points;  Ffix  = St.ConnectivityList;   % tibia (fixed in STL coords)

%% ---------------- MASTER BODY (TIBIA) ----------------
master = buildBodyStruct(Ffix, Vfix);

% Precompute fields used by your hybrid/ROI contact function(s)
master.faceNormals = faceNormalsFromFV(master.F, master.V);
master.kdtreeFaces = KDTreeSearcher(master.centroids);
master.bbox        = bboxFromV(master.V);

%% ---------------- FLEXION-ONLY ALIGNMENT ----------------
flex0 = getFirstFrameFlexion(tr);
Rcorr = rotxd(-flex0);
fprintf('Flexion offset removed: %.4f deg\n', flex0);

%% ---------------- CONSTANT REGISTRATION (FRAME 1) ----------------
% This is the alignment fix: we map STL femur coords to tracking femur coords (with flexion correction).
% Same idea as your "working" animation script.
tTf1 = inv(fTt(:,:,1));                 % femur-from-tibia? (your convention)
tTf1(1:3,1:3) = Rcorr * tTf1(1:3,1:3);  % flexion correction applied to rotation

Tstl = stlFrameFromPCA(Vmov0);          % STL principal axes frame
Treg = tTf1 / Tstl;                     % constant STL->tracking registration

fprintf('Computed constant STL ↔ tracking registration (frame 1).\n');

if debugFrame1Align
    figure('Color','w'); hold on; axis equal; grid on; view(3);
    title('DEBUG: Frame 1 alignment (Treg + flexion correction)');
    xlabel('X'); ylabel('Y'); zlabel('Z');

    patch('Faces',master.F,'Vertices',master.V, ...
        'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',0.25,'EdgeColor','none');

    Vm_dbg = applyTformToVertices(Treg * tTf1, Vmov0);
    patch('Faces',Fmov0,'Vertices',Vm_dbg, ...
        'FaceColor',[0 0.7 0], 'FaceAlpha',0.5,'EdgeColor','none');

    camlight; lighting gouraud
end

%% ---------------- INITIAL MESH ----------------
Flevel = Fmov0;
Vlevel = Vmov0;

%% ---------------- LOGS ----------------
Amean      = nan(1,maxLevels);
pctDiff    = nan(1,maxLevels);
Nfaces     = nan(1,maxLevels);
belowCount = 0;

%% ================= MAIN LOOP =================
for lvl = 1:maxLevels
    fprintf('\n--- LEVEL %d ---\n', lvl);

    Nfaces(lvl) = size(Flevel,1);
    areas = zeros(nFrames,1);

    unionMask = false(size(Flevel,1),1);

    % Precompute adjacency ONCE per level (used for ring growth)
    adj = buildFaceAdjacency_edges(Flevel);

    % ---- CONTACT OVER ALL FRAMES ----
    for k = 1:nFrames
        tTf = inv(fTt(:,:,k));
        tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);

        % Apply the SAME alignment pipeline as the working animation:
        % STL femur -> (Treg * tTf) into tibia/master frame
        T  = Treg * tTf;
        Vm = applyTformToVertices(T, Vlevel);

        slave = buildBodyStruct(Flevel, Vm);
        res   = callContact(master, slave, tol, preferROI);

        areas(k) = res.contactArea;

        if isfield(res,'contactMask') && numel(res.contactMask)==size(Flevel,1)
            unionMask = unionMask | logical(res.contactMask(:));
        end
    end

    Amean(lvl) = mean(areas);

    if lvl > 1
        denom = max(0.5*(Amean(lvl)+Amean(lvl-1)), 1e-6);
        pctDiff(lvl) = 100 * abs(Amean(lvl)-Amean(lvl-1)) / denom;
    end

    fprintf('faces=%d | mean(A)=%.4f | %%diff=%.4f\n', ...
        Nfaces(lvl), Amean(lvl), pctDiff(lvl));

    % ---- OPTIONAL ANIMATE THIS LEVEL ----
    if animateEveryLevel
        animateLevel(lvl, master, Flevel, Vlevel, fTt, Rcorr, Treg, tol, framePause, preferROI);
    end

    % ---- CONVERGENCE CHECK ----
    if lvl >= minLevels && ~isnan(pctDiff(lvl)) && pctDiff(lvl) < convTolPct
        belowCount = belowCount + 1;
    else
        belowCount = 0;
    end

    if lvl >= minLevels && belowCount >= convConsec
        fprintf('>>> CONVERGED at level %d <<<\n', lvl);
        break
    end

    % ---- REFINE UNION CONTACT REGION ----
    if ~any(unionMask)
        warning('No contact — stopping.');
        break
    end

    mask = unionMask;
    for r = 1:expandRings
        mask = growMaskOneRing_fast(adj, mask);
    end

    [Flevel, Vlevel] = refineTrianglesMidpoint_cached_simple(Flevel, Vlevel, mask);
end

convergedLevel = lvl;

%% ---- Animate final only (fast workflow) ----
if ~animateEveryLevel && animateFinalOnly
    animateLevel(convergedLevel, master, Flevel, Vlevel, fTt, Rcorr, Treg, tol, framePause, preferROI);
end

%% ---------------- EXPORT CONVERGED STL ----------------
outDir = fullfile(projectRoot,'results_convergedSTL');
if ~exist(outDir,'dir'), mkdir(outDir); end

outStl = fullfile(outDir, sprintf('femur_converged_level_%02d.stl', convergedLevel));
stlwrite(outStl, Flevel, Vlevel);   % exports in *LOCAL* femur STL coordinates (as requested)

save(fullfile(outDir,'convergence_meta.mat'), ...
    'Amean','pctDiff','Nfaces','convergedLevel','tol','expandRings');

fprintf('\nSaved converged STL:\n  %s\n', outStl);

%% ---------------- PLOT CONVERGENCE ----------------
figure('Color','w'); grid on; hold on;
yyaxis left;  plot(Amean,'-o','LineWidth',1.5); ylabel('Mean contact area (mm^2)');
yyaxis right; plot(pctDiff,'-s','LineWidth',1.5); ylabel('% diff vs previous level');
xlabel('Refinement level');
title('Dynamic local mesh convergence (mean area across frames)');
set(gca,'XLim',[1 maxLevels]);

%% =======================================================================
%%                               HELPERS
%% =======================================================================

function res = callContact(master, slave, tol, preferROI)
    if preferROI && exist('computeContactArea_STS_hybrid_ROI','file')==2
        res = computeContactArea_STS_hybrid_ROI(master, slave, tol);
    else
        res = computeContactArea_STS_hybrid(master, slave, tol);
    end
end

function flex0 = getFirstFrameFlexion(tr)
    TF = tr.Kinematics.tibiofemoral;
    if istable(TF)
        flex0 = TF{1,'flexion'};
    else
        flex0 = TF(1,1);
    end
end

function body = buildBodyStruct(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    body.F = F;
    body.V = V;
    body.centroids = (v1+v2+v3)/3;
    body.areas = 0.5*vecnorm(cross(v2-v1,v3-v1,2),2,2);
end

function Vout = applyTformToVertices(T,V)
    Vh = [V ones(size(V,1),1)];
    Vt = (T*Vh')';
    Vout = Vt(:,1:3);
end

function R = rotxd(deg)
    a = deg2rad(deg);
    R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end

function Tstl = stlFrameFromPCA(V)
    c = mean(V,1);
    [~,~,Vh] = svd(V - c, 'econ');
    R = Vh;
    if det(R) < 0
        R(:,3) = -R(:,3);
    end
    Tstl = eye(4);
    Tstl(1:3,1:3) = R;
    Tstl(1:3,4)   = c(:);
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
end

function N = faceNormalsFromFV(F, V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    N = cross(v2 - v1, v3 - v1, 2);
    mag = sqrt(sum(N.^2,2));
    mag(mag==0) = 1;
    N = N ./ mag;
end

%% ---------- FAST mask growth via FACE adjacency (shared edges) ----------

function adj = buildFaceAdjacency_edges(F)
% adj{i} = faces that share an edge with face i (edge-based adjacency)
    nF = size(F,1);
    adj = cell(nF,1);

    E = [F(:,[1 2]);
         F(:,[2 3]);
         F(:,[3 1])];
    faceID = repelem((1:nF)',3);

    E = sort(E,2);
    [~,~,ic] = unique(E,'rows');

    for k = 1:max(ic)
        faces = faceID(ic==k);
        if numel(faces) > 1
            for f = faces'
                adj{f} = unique([adj{f}; faces(faces~=f)]);
            end
        end
    end
end

function mask2 = growMaskOneRing_fast(adj, mask)
    mask = logical(mask(:));
    mask2 = mask;
    idx = find(mask);
    for ii = 1:numel(idx)
        f = idx(ii);
        mask2(adj{f}) = true;
    end
end

%% ---------- LOCAL midpoint refinement (cached edge midpoints) ----------

function [Fnew,Vnew] = refineTrianglesMidpoint_cached_simple(F,V,mask)
    mask = logical(mask(:));
    Vnew = V;
    edgeMid = containers.Map('KeyType','char','ValueType','int32');

    % Preallocate: each refined tri becomes 4 => +3 faces per refined face
    Fnew = zeros(size(F,1) + 3*nnz(mask), 3);
    o = 0;

    for i = 1:size(F,1)
        a = F(i,1); b = F(i,2); c = F(i,3);
        if mask(i)
            ab = mid(a,b); bc = mid(b,c); ca = mid(c,a);

            o=o+1; Fnew(o,:) = [a  ab ca];
            o=o+1; Fnew(o,:) = [b  bc ab];
            o=o+1; Fnew(o,:) = [c  ca bc];
            o=o+1; Fnew(o,:) = [ab bc ca];
        else
            o=o+1; Fnew(o,:) = [a b c];
        end
    end
    Fnew = Fnew(1:o,:);

    function m = mid(i1,i2)
        k = sprintf('%d_%d', min(i1,i2), max(i1,i2));
        if isKey(edgeMid,k)
            m = edgeMid(k);
        else
            Vnew(end+1,:) = 0.5*(Vnew(i1,:) + Vnew(i2,:));
            m = int32(size(Vnew,1));
            edgeMid(k) = m;
        end
    end
end

%% ---------- Animation (optional) ----------

function animateLevel(lvl, master, Fmov, Vmov, fTt, Rcorr, Treg, tol, pauseT, preferROI)
    figId = 700;
    figure(figId); clf; hold on; axis equal; grid on; view(3);
    title(sprintf('Level %d | faces=%d', lvl, size(Fmov,1)));
    xlabel('X'); ylabel('Y'); zlabel('Z');

    patch('Faces',master.F,'Vertices',master.V,'FaceAlpha',0.2,'EdgeColor','none', ...
        'FaceColor',[0.7 0.7 0.7]);

    h = patch('Faces',Fmov,'Vertices',Vmov, ...
        'FaceColor','flat','EdgeColor','none','FaceAlpha',0.8);

    colormap([0.7 0.7 0.7; 1 0 1]); caxis([1 2]);

    for k = 1:size(fTt,3)
        tTf = inv(fTt(:,:,k));
        tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);

        T  = Treg * tTf;
        Vk = applyTformToVertices(T, Vmov);

        set(h,'Vertices',Vk);

        res = callContact(master, buildBodyStruct(Fmov, Vk), tol, preferROI);

        C = ones(size(Fmov,1),1);
        if isfield(res,'contactMask') && numel(res.contactMask)==size(Fmov,1)
            C(logical(res.contactMask)) = 2;
        end
        set(h,'FaceVertexCData',C);

        drawnow;
        pause(pauseT);
    end
end
