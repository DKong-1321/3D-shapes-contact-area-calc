%% cylinderDynamicConvergence_meanArea_showEveryLevel.m
% -------------------------------------------------------------------------
% Dynamic mesh convergence for meshed CYLINDERS
%
% Convergence metric:
%   - Use OVERALL MEAN contact area per refinement level
%   - % difference computed using symmetric, stable formulation
%   - Stop when mean % diff < 5% for 3 consecutive levels
%
% Visualisation:
%   - Animate ALL frames at EVERY refinement level
%
% Alignment:
%   - FLEXION ONLY correction from first kinematics row
%   - NO anterior modification
%   - Applied correctly to femur motion (tTf)
% -------------------------------------------------------------------------

clear; clc; close all;

% PATHS 
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;
addpath(genpath(fullfile(projectRoot,'src')));

dataDir  = fullfile(projectRoot,'data');
modelDir = fullfile(projectRoot,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');

assert(exist(matPath,'file')==2);
assert(exist(femurSTL,'file')==2);
assert(exist(tibiaSTL,'file')==2);

% OPTIONS
trialIndex = 1;
tol        = 0.20;     % mm
maxLevels  = 25;
minLevels  = 3;

convTolPct = 5;        
convConsec = 3;

expandRings = 1;       % expand union contact region before refinement
framePause  = 0.02;

% LOAD DATA
S  = load(matPath);
tr = S.trajectories(trialIndex);

fTt = tr.Transform.fTt;
nFrames = size(fTt,3);

% LOAD STLs 
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm0 = Sf.ConnectivityList;
Vt  = St.Points;  Ft  = St.ConnectivityList;

% MASTER BODY (tibia fixed)
master = buildBodyStruct(Ft, Vt);

% FLEXION-ONLY ALIGNMENT FIX
flex0 = getFirstFrameFlexion(tr);
fprintf('Initial flexion offset = %.4f deg\n', flex0);

Rcorr = rotxd(-flex0);   % remove offset

% INITIAL MESH
Flevel = Fm0;
Vlevel = Vm0;

% LOGS
Amean   = nan(1,maxLevels);
pctDiff = nan(1,maxLevels);
Nfaces  = nan(1,maxLevels);

belowCount = 0;

fprintf('\n=== MEAN-AREA CONVERGENCE ===\n');

% MAIN REFINEMENT LOOP
for lvl = 1:maxLevels

    fprintf('\n--- LEVEL %d ---\n', lvl);
    Nfaces(lvl) = size(Flevel,1);

    areas = zeros(nFrames,1);
    unionMask = false(size(Flevel,1),1);

    % CONTACT OVER ALL FRAMES
    for k = 1:nFrames
        tTf = inv(fTt(:,:,k));
        tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);

        Vm_k = applyTformToVertices(tTf, Vlevel);
        slave = buildBodyStruct(Flevel, Vm_k);

        res = computeContactArea_STS_hybrid(master, slave, tol);
        areas(k) = res.contactArea;

        if isfield(res,'contactMask') && numel(res.contactMask)==size(Flevel,1)
            unionMask = unionMask | res.contactMask;
        end
    end

    % MEAN AREA
    Amean(lvl) = mean(areas);

    % DIFFERENCE (STABLE)
    if lvl > 1
        denom = max(0.5*(Amean(lvl)+Amean(lvl-1)), 1e-6);
        pctDiff(lvl) = 100 * abs(Amean(lvl)-Amean(lvl-1)) / denom;
    end

    fprintf('faces=%d | mean(A)=%.4f | %%diff=%.4f\n', ...
        Nfaces(lvl), Amean(lvl), pctDiff(lvl));

    % ANIMATE THIS LEVEL
    animateLevel(lvl, master, Flevel, Vlevel, fTt, Rcorr, tol, framePause);

    % CONVERGENCE CHECK 
    if lvl >= minLevels && ~isnan(pctDiff(lvl)) && pctDiff(lvl) < convTolPct
        belowCount = belowCount + 1;
    else
        belowCount = 0;
    end

    if lvl >= minLevels && belowCount >= convConsec
        fprintf('>>> CONVERGED at level %d <<<\n', lvl);
        break
    end

    % REFINE UNION CONTACT REGION 
    if ~any(unionMask)
        warning('No contact at level %d — stopping.', lvl);
        break
    end

    mask = unionMask;
    for r = 1:expandRings
        mask = expandOneRing(Flevel, mask);
    end

    [Flevel, Vlevel] = refineTrianglesMidpoint_cached_simple(Flevel, Vlevel, mask);
end

% CONVERGENCE PLOT
figure('Color','w'); grid on; hold on;
plot(1:lvl, Amean(1:lvl), '-o');
xlabel('Refinement level');
ylabel('Mean contact area (mm^2)');
title('Mean contact area convergence');

figure('Color','w'); grid on; hold on;
plot(1:lvl, pctDiff(1:lvl), '-o');
yline(convTolPct,'--');
xlabel('Refinement level');
ylabel('% difference');
title('Convergence criterion');

fprintf('\n=== DONE ===\n');

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

function Vout = applyTformToVertices(T, V)
    Vh = [V ones(size(V,1),1)];
    Vt = (T * Vh')';
    Vout = Vt(:,1:3);
end

function R = rotxd(deg)
    a = deg2rad(deg);
    R = [1 0 0;
         0 cos(a) -sin(a);
         0 sin(a)  cos(a)];
end

function mask2 = expandOneRing(F, mask)
    mask2 = mask;
    nV = max(F(:));
    v2f = cell(nV,1);
    for fi = 1:size(F,1)
        for j = 1:3
            v2f{F(fi,j)}(end+1) = fi;
        end
    end
    for fi = find(mask)'
        for v = F(fi,:)
            mask2(v2f{v}) = true;
        end
    end
end

function [Fnew, Vnew] = refineTrianglesMidpoint_cached_simple(F, V, mask)
    Vnew = V;
    nF = size(F,1);
    Fnew = zeros(nF + 3*nnz(mask),3);
    edgeMid = containers.Map('KeyType','char','ValueType','int32');

    o = 0;
    for fi = 1:nF
        a = F(fi,1); b = F(fi,2); c = F(fi,3);
        if mask(fi)
            ab = mid(a,b); bc = mid(b,c); ca = mid(c,a);
            o=o+1; Fnew(o,:)=[a ab ca];
            o=o+1; Fnew(o,:)=[b bc ab];
            o=o+1; Fnew(o,:)=[c ca bc];
            o=o+1; Fnew(o,:)=[ab bc ca];
        else
            o=o+1; Fnew(o,:)=[a b c];
        end
    end
    Fnew = Fnew(1:o,:);

    function m = mid(i,j)
        key = sprintf('%d_%d',min(i,j),max(i,j));
        if isKey(edgeMid,key)
            m = edgeMid(key);
        else
            Vnew(end+1,:) = 0.5*(Vnew(i,:) + Vnew(j,:));
            m = size(Vnew,1);
            edgeMid(key) = m;
        end
    end
end

function animateLevel(lvl, master, Fmov, Vmov, fTt, Rcorr, tol, pauseT)
    figure(800+lvl); clf; hold on; axis equal; grid on; view(3);
    title(sprintf('Level %d | faces=%d', lvl, size(Fmov,1)));
    patch('Faces',master.F,'Vertices',master.V,'FaceAlpha',0.2,'EdgeColor','none');

    h = patch('Faces',Fmov,'Vertices',Vmov,'FaceColor','flat','EdgeColor','none');
    colormap([0.7 0.7 0.7; 1 0 1]); caxis([1 2]);

    for k = 1:size(fTt,3)
        tTf = inv(fTt(:,:,k));
        tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);
        Vk = applyTformToVertices(tTf, Vmov);
        set(h,'Vertices',Vk);

        res = computeContactArea_STS_hybrid(master, buildBodyStruct(Fmov,Vk), tol);
        C = ones(size(Fmov,1),1);
        if isfield(res,'contactMask') && numel(res.contactMask)==size(Fmov,1)
            C(res.contactMask) = 2;
        end
        set(h,'FaceVertexCData',C);

        drawnow;
        pause(pauseT);
    end
end
