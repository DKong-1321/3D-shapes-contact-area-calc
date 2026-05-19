%% cylinderConvergence_meanArea_exportSTL.m
% -------------------------------------------------------------------------
% Dynamic mesh convergence for meshed CYLINDERS
%
% - Uses OVERALL MEAN contact area per refinement level
% - Stable symmetric % difference metric
% - Converges when mean % diff < 5% for 3 consecutive levels
% - Animates ALL frames at EVERY refinement level
% - Exports converged FEMUR mesh as STL (local coordinates)
% - Flexion-only alignment (NO anterior)
% -------------------------------------------------------------------------

clear; clc; close all;

%% ---------------- PATHS ----------------
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

%% ---------------- OPTIONS ----------------
trialIndex = 1;
tol        = 0.20;      % mm
maxLevels  = 25;
minLevels  = 3;
convTolPct = 5;
convConsec = 3;
expandRings = 1;
framePause = 0.02;

%% ---------------- LOAD DATA ----------------
S  = load(matPath);
tr = S.trajectories(trialIndex);
fTt = tr.Transform.fTt;
nFrames = size(fTt,3);

%% ---------------- LOAD STLs ----------------
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vmov0 = Sf.Points;  Fmov0 = Sf.ConnectivityList;
Vfix  = St.Points;  Ffix  = St.ConnectivityList;

%% ---------------- MASTER BODY ----------------
master = buildBodyStruct(Ffix, Vfix);

%% ---------------- FLEXION-ONLY ALIGNMENT ----------------
flex0 = getFirstFrameFlexion(tr);
Rcorr = rotxd(-flex0);
fprintf('Flexion offset removed: %.4f deg\n', flex0);

%% ---------------- INITIAL MESH ----------------
Flevel = Fmov0;
Vlevel = Vmov0;

%% ---------------- LOGS ----------------
Amean   = nan(1,maxLevels);
pctDiff = nan(1,maxLevels);
Nfaces  = nan(1,maxLevels);
belowCount = 0;

%% ================= MAIN LOOP =================
for lvl = 1:maxLevels
    fprintf('\n--- LEVEL %d ---\n', lvl);

    Nfaces(lvl) = size(Flevel,1);
    areas = zeros(nFrames,1);
    unionMask = false(size(Flevel,1),1);

    %% ---- CONTACT OVER ALL FRAMES ----
    for k = 1:nFrames
        tTf = inv(fTt(:,:,k));
        tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);

        Vm = applyTformToVertices(tTf, Vlevel);
        slave = buildBodyStruct(Flevel, Vm);

        res = computeContactArea_STS_hybrid(master, slave, tol);
        areas(k) = res.contactArea;

        if isfield(res,'contactMask') && numel(res.contactMask)==size(Flevel,1)
            unionMask = unionMask | res.contactMask;
        end
    end

    Amean(lvl) = mean(areas);

    if lvl > 1
        denom = max(0.5*(Amean(lvl)+Amean(lvl-1)), 1e-6);
        pctDiff(lvl) = 100 * abs(Amean(lvl)-Amean(lvl-1)) / denom;
    end

    fprintf('faces=%d | mean(A)=%.4f | %%diff=%.4f\n', ...
        Nfaces(lvl), Amean(lvl), pctDiff(lvl));

    %% ---- ANIMATE THIS LEVEL ----
    animateLevel(lvl, master, Flevel, Vlevel, fTt, Rcorr, tol, framePause);

    %% ---- CONVERGENCE CHECK ----
    if lvl >= minLevels && ~isnan(pctDiff(lvl)) && pctDiff(lvl) < convTolPct
        belowCount = belowCount + 1;
    else
        belowCount = 0;
    end

    if lvl >= minLevels && belowCount >= convConsec
        fprintf('>>> CONVERGED at level %d <<<\n', lvl);
        break
    end

    %% ---- REFINE UNION CONTACT REGION ----
    if ~any(unionMask)
        warning('No contact — stopping.');
        break
    end

    mask = unionMask;
    for r = 1:expandRings
        mask = expandOneRing(Flevel, mask);
    end

    [Flevel, Vlevel] = refineTrianglesMidpoint_cached_simple(Flevel, Vlevel, mask);
end

convergedLevel = lvl;

%% ---------------- EXPORT CONVERGED STL ----------------
outDir = fullfile(projectRoot,'results_convergedSTL');
if ~exist(outDir,'dir'), mkdir(outDir); end

outStl = fullfile(outDir, sprintf('femur_converged_level_%02d.stl', convergedLevel));
stlwrite(outStl, Flevel, Vlevel);

save(fullfile(outDir,'convergence_meta.mat'), ...
    'Amean','pctDiff','Nfaces','convergedLevel','tol');

fprintf('\nSaved converged STL:\n  %s\n', outStl);

%% =======================================================================
%%                               HELPERS
%% =======================================================================

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

function mask2 = expandOneRing(F,mask)
    mask2 = mask;
    v2f = cell(max(F(:)),1);
    for i=1:size(F,1)
        for j=1:3
            v2f{F(i,j)}(end+1)=i;
        end
    end
    for i=find(mask)'
        for v=F(i,:)
            mask2(v2f{v})=true;
        end
    end
end

function [Fnew,Vnew] = refineTrianglesMidpoint_cached_simple(F,V,mask)
    Vnew = V;
    edgeMid = containers.Map('KeyType','char','ValueType','int32');
    Fnew = zeros(size(F,1)+3*nnz(mask),3);
    o=0;
    for i=1:size(F,1)
        a=F(i,1); b=F(i,2); c=F(i,3);
        if mask(i)
            ab=mid(a,b); bc=mid(b,c); ca=mid(c,a);
            o=o+1; Fnew(o,:)=[a ab ca];
            o=o+1; Fnew(o,:)=[b bc ab];
            o=o+1; Fnew(o,:)=[c ca bc];
            o=o+1; Fnew(o,:)=[ab bc ca];
        else
            o=o+1; Fnew(o,:)=[a b c];
        end
    end
    Fnew=Fnew(1:o,:);
    function m=mid(i1,i2)
        k=sprintf('%d_%d',min(i1,i2),max(i1,i2));
        if isKey(edgeMid,k)
            m=edgeMid(k);
        else
            Vnew(end+1,:)=0.5*(Vnew(i1,:)+Vnew(i2,:));
            m=size(Vnew,1);
            edgeMid(k)=m;
        end
    end
end

function animateLevel(lvl, master, Fmov, Vmov, fTt, Rcorr, tol, pauseT)
    figure(700+lvl); clf; hold on; axis equal; grid on; view(3);
    title(sprintf('Level %d | faces=%d',lvl,size(Fmov,1)));
    patch('Faces',master.F,'Vertices',master.V,'FaceAlpha',0.2,'EdgeColor','none');
    h=patch('Faces',Fmov,'Vertices',Vmov,'FaceColor','flat','EdgeColor','none');
    colormap([0.7 0.7 0.7;1 0 1]); caxis([1 2]);
    for k=1:size(fTt,3)
        tTf=inv(fTt(:,:,k));
        tTf(1:3,1:3)=Rcorr*tTf(1:3,1:3);
        Vk=applyTformToVertices(tTf,Vmov);
        set(h,'Vertices',Vk);
        res=computeContactArea_STS_hybrid(master,buildBodyStruct(Fmov,Vk),tol);
        C=ones(size(Fmov,1),1);
        if isfield(res,'contactMask') && numel(res.contactMask)==size(Fmov,1)
            C(res.contactMask)=2;
        end
        set(h,'FaceVertexCData',C);
        drawnow; pause(pauseT);
    end
end
