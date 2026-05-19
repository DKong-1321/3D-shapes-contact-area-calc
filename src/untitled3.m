% sphereCube_localRefine_contactTrisOnly_autoStop.m
% ------------------------------------------------------------
% Sphere–Cube local surface refinement using ONLY the triangles
% used in contact on the previous step, refining until contact
% area change is < 5% (for N consecutive steps).
%
% Outputs:
%  1) Coarse contact (mesh edges visible)
%  2) Step-before-final contact
%  3) Final contact
%  4) Final "zoned refinement" plot (faces coloured by refinement level)
%
% REQUIREMENTS:
%  - computeContactArea_STS
%  - buildBodyStruct
%  - loadStlMesh OR loadAnyStl
% ------------------------------------------------------------

clear; clc; close all;

%% --- RENDERING (Mac-friendly) ---
% Renderer preference warnings can be ignored; keep it simple:
opengl('save','software');  % forces software OpenGL (avoids WebGL crashes)

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')));
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Sphere–Cube LOCAL CONTACT REFINEMENT (auto-stop) ===\n');
fprintf('Project root: %s\n',projectRoot);
fprintf('Model dir:    %s\n\n',modelDir);

%% ===== FILES =====
cubeStl   = fullfile(modelDir,'50mm cube 4 iterations.stl');      % MASTER
sphereStl = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');   % MOVING

assert(exist(cubeStl,'file')==2,'Cube STL missing: %s', cubeStl);
assert(exist(sphereStl,'file')==2,'Sphere STL missing: %s', sphereStl);

%% ===== CONTACT SETTINGS =====
tol = 0.20;                % mm
penetration = 0.5*tol;      % guaranteed contact

opts = struct();
opts.tol               = tol;
opts.roiExpandFactor   = 1.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 20;
opts.sampleMode        = 'centroid';
opts.sampleThreshold   = 0.5;
opts.useCentroidPrefilter = true;
opts.prefilterFactor      = 6.0;
opts.returnMask = true;    % try to request contact triangle mask

%% ===== LOCAL REFINEMENT CONTROL =====
HmaxCoarse  = 6.0;          % mm initial sphere mesh (generateMesh surface)

% ---- AUTO STOP SETTINGS ----
maxSteps        = 25;       % safety cap
targetPct       = 5.0;      % stop when % change in area < 5%
needConsecutive = 2;        % require this many consecutive steps under target
minValidArea    = 1e-12;    % ignore 0/near-0 areas in convergence check

%% ===== LOAD MASTER CUBE =====
if ~exist('buildBodyStruct','file'), error('buildBodyStruct not found on path.'); end
if ~exist('computeContactArea_STS','file'), error('computeContactArea_STS not found on path.'); end

[Fcube,Vcube] = loadMeshAny(cubeStl);
cubeBody = buildBodyStruct(Fcube,Vcube);
cubeBody.bbox  = bboxFromV(Vcube);
cubeBody.kdtree = KDTreeSearcher(cubeBody.centroids);

%% ===== INITIAL COARSE SPHERE SURFACE =====
[Fmov,Vmov] = surfaceFromStl_generateMesh(sphereStl,HmaxCoarse);

% refinement level per face (0 = coarse)
faceLvl = zeros(size(Fmov,1),1,'uint16');

% Place sphere on cube top
[VcubeP,VmovP] = placeSphereOnCubeTop(Vcube,Vmov,penetration);
cubeBody.V = VcubeP;

%% ===== ITERATIVE CONTACT → LOCAL REFINE (AUTO STOP) =====
snap = struct('F',{},'V',{},'mask',{},'A',{},'faceLvl',{},'dPct',{});

belowCount = 0;
Aprev = NaN;

for s = 1:maxSteps
    movBody = buildBodyStruct(Fmov, VmovP);
    movBody.bbox = bboxFromV(VmovP);

    % --- Run contact ---
    res = computeContactArea_STS(cubeBody, movBody, tol, opts);

    A = NaN; mask = [];
    if isstruct(res)
        if isfield(res,'contactArea'), A = res.contactArea; end
        if isfield(res,'contactMask'), mask = res.contactMask; end
        if isempty(mask) && isfield(res,'maskSlave'), mask = res.maskSlave; end
        if isempty(mask) && isfield(res,'mask'), mask = res.mask; end
    else
        A = res;
    end

    % --- Fallback mask (only if computeContactArea_STS didn't return one) ---
    if isempty(mask)
        mask = fallbackMaskByDistance(cubeBody, movBody, tol, 1.0);
        warning('Fallback mask used (distance-based).');
    end

    % --- % change vs previous valid area ---
    dPct = NaN;
    if ~isnan(Aprev) && ~isnan(A) && Aprev > minValidArea
        dPct = 100 * abs(A - Aprev) / Aprev;
    end

    fprintf('Step %d/%d | faces=%d | contact=%d | A=%.6f mm^2 | dA=%.2f%%\n', ...
        s, maxSteps, size(Fmov,1), nnz(mask), A, dPct);

    % --- store snapshot ---
    snap(s).F = Fmov;
    snap(s).V = VmovP;
    snap(s).mask = mask;
    snap(s).A = A;
    snap(s).faceLvl = faceLvl;
    snap(s).dPct = dPct;

    % --- stop if no contact ---
    if ~any(mask)
        warning('No contact triangles detected at step %d. Stopping refinement.', s);
        break
    end

    % --- stopping rule: below targetPct for needConsecutive steps ---
    if ~isnan(dPct) && dPct < targetPct
        belowCount = belowCount + 1;
    else
        belowCount = 0;
    end

    % update Aprev only when A is valid and > 0
    if ~isnan(A) && A > minValidArea
        Aprev = A;
    end

    if belowCount >= needConsecutive
        fprintf('STOP: dA < %.2f%% for %d consecutive steps.\n', targetPct, needConsecutive);
        break
    end

    % --- refine ONLY contact triangles ---
    [Fmov, VmovP, faceLvl] = refineTrianglesMidpoint_cached(Fmov, VmovP, mask, faceLvl);
end

% Trim snap if loop ended early
snap = snap(1:min(numel(snap), s));

%% ===== VISUALISATIONS =====
plotStep(snap, 1, cubeBody, tol, 'COARSE');

nSnap = numel(snap);
if nSnap >= 2
    plotStep(snap, max(1, nSnap-1), cubeBody, tol, 'STEP BEFORE FINAL');
end

plotStep(snap, nSnap, cubeBody, tol, 'FINAL CONTACT');

plotFinalZoning(snap(nSnap), cubeBody, tol);

fprintf('\n=== DONE ===\n');

%% ============================================================
%% ====================== FUNCTIONS ============================
%% ============================================================

function [F,V] = loadMeshAny(path)
    if exist('loadStlMesh','file')
        [F,V] = loadStlMesh(path);
    elseif exist('loadAnyStl','file')
        [F,V] = loadAnyStl(path);
    else
        error('No STL loader found (need loadStlMesh or loadAnyStl).');
    end
end

function [Fsurf,Vsurf] = surfaceFromStl_generateMesh(stlPath,Hmax)
    model = createpde();
    importGeometry(model, stlPath);
    try
        msh = generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');
    catch
        msh = generateMesh(model,'Hmax',Hmax);
    end
    V = msh.Nodes.';
    E = msh.Elements;
    if size(E,1) < 4
        % already triangles
        Fsurf = E.';
        Vsurf = V;
        return
    end
    T = E(1:4,:).';  % take corner nodes
    TR = triangulation(T,V);
    [Fsurf,Vsurf] = freeBoundary(TR);
end

function [VcubeP,VsphP] = placeSphereOnCubeTop(Vcube,Vsph,penetration)
    VcubeP = Vcube;

    cC = mean(Vcube(:,1:2),1);
    cS = mean(Vsph(:,1:2),1);

    VsphP = Vsph;
    VsphP(:,1:2) = VsphP(:,1:2) + (cC - cS);

    zTop = max(VcubeP(:,3));
    zBot = min(VsphP(:,3));
    VsphP(:,3) = VsphP(:,3) + ((zTop - penetration) - zBot);
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1));
          min(V(:,2)) max(V(:,2));
          min(V(:,3)) max(V(:,3))];
end

function mask = fallbackMaskByDistance(master,mov,tol,factor)
    if ~isfield(master,'kdtree') || isempty(master.kdtree)
        master.kdtree = KDTreeSearcher(master.centroids);
    end
    [~,d] = knnsearch(master.kdtree, mov.centroids, 'K', 1);
    mask = d <= factor * tol;
end

function [F2,V2,lvl2] = refineTrianglesMidpoint_cached(F,V,mask,lvl)
% Midpoint subdivision ONLY on faces where mask==true.
% Tracks refinement level per face:
% - kept faces keep their level
% - refined faces produce 4 children at parentLevel+1

    idxR = find(mask);
    idxK = find(~mask);

    Fk = F(idxK,:);  lvlK = lvl(idxK);
    Fr = F(idxR,:);  lvlR = lvl(idxR);

    V2 = V;

    edgeMap = containers.Map('KeyType','char','ValueType','int32');

    Fnew   = zeros(4*numel(idxR),3);
    lvlNew = zeros(4*numel(idxR),1,'uint16');

    o = 0;
    for i = 1:size(Fr,1)
        a = Fr(i,1); b = Fr(i,2); c = Fr(i,3);

        iab = mid(a,b);
        ibc = mid(b,c);
        ica = mid(c,a);

        L = uint16(lvlR(i) + 1);

        o=o+1; Fnew(o,:) = [a,   iab, ica]; lvlNew(o) = L;
        o=o+1; Fnew(o,:) = [iab, b,   ibc]; lvlNew(o) = L;
        o=o+1; Fnew(o,:) = [ica, ibc, c  ]; lvlNew(o) = L;
        o=o+1; Fnew(o,:) = [iab, ibc, ica]; lvlNew(o) = L;
    end

    F2   = [Fk; Fnew(1:o,:)];
    lvl2 = [lvlK; lvlNew(1:o)];

    [V2,F2,lvl2] = removeUnusedVertices_keepFaceLevels(V2,F2,lvl2);

    function m = mid(i1,i2)
        k = sprintf('%d_%d', min(i1,i2), max(i1,i2));
        if isKey(edgeMap,k)
            m = edgeMap(k);
            return
        end
        V2(end+1,:) = 0.5*(V2(i1,:) + V2(i2,:));
        m = size(V2,1);
        edgeMap(k) = m;
    end
end

function [V2,F2,lvl2] = removeUnusedVertices_keepFaceLevels(V,F,lvl)
    used = false(size(V,1),1);
    used(F(:)) = true;

    map = zeros(size(V,1),1);
    map(used) = 1:nnz(used);

    V2 = V(used,:);
    F2 = map(F);

    lvl2 = lvl; % face list unchanged by vertex compaction
end

function plotStep(snap,idx,master,tol,label)
    F = snap(idx).F;
    V = snap(idx).V;
    mask = snap(idx).mask;
    A = snap(idx).A;

    figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | tol=%.3f mm | A=%.6f mm^2', label, tol, A), 'Interpreter','none');

    % Master cube (faint)
    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    % Sphere mesh (edges visible)
    patch('Faces', F, 'Vertices', V, ...
        'FaceColor', [0 .7 0], 'FaceAlpha', 0.08, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, ...
        'LineWidth', 0.25);

    % Contact triangles highlighted
    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], 'FaceAlpha', 0.90, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 0.25, ...
            'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end

function plotFinalZoning(snap,master,tol)
    F = snap.F;
    V = snap.V;
    lvl = double(snap.faceLvl);
    mask = snap.mask;
    A = snap.A;

    figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('FINAL ZONED LOCAL REFINEMENT | tol=%.3f mm | A=%.6f mm^2', tol, A), 'Interpreter','none');

    % Master cube faint
    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    % Sphere coloured by refinement level (shows rings/bands)
    patch('Faces', F, 'Vertices', V, ...
        'FaceVertexCData', lvl, ...
        'FaceColor', 'flat', ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', [0 0 0], ...
        'EdgeAlpha', 0.22, ...
        'LineWidth', 0.25);

    cb = colorbar;
    cb.Label.String = 'Refinement level (0=coarse, higher=finer)';

    % Overlay final contact region
    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], ...
            'FaceAlpha', 0.95, ...
            'EdgeColor', [0 0 0], ...
            'EdgeAlpha', 0.25, ...
            'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end
