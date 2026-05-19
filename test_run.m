% sphereCube_localRefine_contactTrisOnly.m
% ------------------------------------------------------------
% Sphere–Cube local surface refinement using ONLY the triangles
% that were used to calculate contact area on the previous step.
%
% Workflow (what you asked for):
%  1) Build COARSE moving sphere surface mesh (generateMesh -> freeBoundary)
%  2) Place sphere on top of cube (small penetration)
%  3) Run contact area -> get "contact triangles mask" on the sphere
%  4) Refine ONLY those triangles on the sphere (midpoint subdivision)
%  5) Re-run contact area and repeat
%
% At the end it shows:
%  - COARSE step (contact triangles highlighted)
%  - STEP BEFORE FINAL (highlighted)
%  - FINAL step (highlighted)
%
% Notes:
%  - This assumes computeContactArea_STS can return a contact triangle mask
%    when opts.returnMask=true (field names vary). The script tries common
%    names: contactMask, maskSlave, mask.
%  - If no mask is returned, it falls back to an approximate mask based on
%    distance-to-master centroids (still "area-of-interest" but not exact).
%  - Plotting is face-only + decimated for display to avoid WebGL crashes.
% ------------------------------------------------------------

clear; clc; close all;

% ---- prevent WebGL crashes on big patches (Mac-friendly) ----
set(groot, 'defaultFigureRenderer', 'opengl');
opengl('save','software');  % forces software OpenGL

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Sphere–Cube local refinement (contact triangles only) ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

%% ===== FILES (EDIT THESE EXACTLY) =====
cubeStl   = fullfile(modelDir,'50mm cube 4 iterations.stl');      % MASTER (fixed)
sphereStl = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');   % MOVING (refined locally)

assert(exist(cubeStl,'file')==2,   'Cube STL not found: %s', cubeStl);
assert(exist(sphereStl,'file')==2, 'Sphere STL not found: %s', sphereStl);

%% ===== CONTACT SETTINGS =====
tol = 0.20;                 % mm
penetration = 0.5*tol;      % small overlap to guarantee contact

opts = struct();
opts.tol               = tol;
opts.roiExpandFactor   = 1.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 20;
opts.sampleMode        = 'centroid';
opts.sampleThreshold   = 0.5;

% request mask (if your computeContactArea_STS supports it)
opts.returnMask        = true;

% any existing prefilter your code uses
opts.useCentroidPrefilter = true;
opts.prefilterFactor      = 6.0;

%% ===== LOCAL REFINEMENT CONTROLS =====
nLocalSteps = 4;      % number of "coarse->refine->rerun" cycles (including coarse)
% Example: 4 means: step1 coarse, step2 refined, step3 refined, step4 refined

% Fallback mask logic (ONLY used if computeContactArea_STS doesn't return a mask)
fallbackDistFactor = 1.0;   % triangles within (factor*tol) from cube surface centroids

%% ===== COARSE SURFACE MESH SETTINGS (generateMesh) =====
HmaxCoarse = 6.0;   % mm start for sphere surface generation (coarse)

%% ===== LOAD MASTER CUBE (fixed) =====
if ~exist('buildBodyStruct','file'), error('buildBodyStruct not found on path.'); end
if ~exist('computeContactArea_STS','file'), error('computeContactArea_STS not found on path.'); end

[Fcube, Vcube] = loadMeshAny(cubeStl);
cubeBody = buildBodyStruct(Fcube, Vcube);
cubeBody.bbox   = bboxFromV(Vcube);
cubeBody.kdtree = KDTreeSearcher(cubeBody.centroids);

%% ===== BUILD INITIAL COARSE MOVING SPHERE SURFACE =====
[Fmov, Vmov] = surfaceFromStl_generateMesh(sphereStl, HmaxCoarse);

% Place sphere on top of cube
[VfixP, VmovP] = placeSphereOnCubeTop(Vcube, Vmov, penetration);

% We keep cube fixed (master), and only move sphere vertices
cubeBodyPlaced = cubeBody;
cubeBodyPlaced.V = VfixP;                % same as Vcube here
cubeBodyPlaced.bbox = bboxFromV(VfixP);  % same

%% ===== ITERATE: contact -> refine only contact tris -> contact -> ... =====
snap = struct('F',{},'V',{},'mask',{},'A',{});

for s = 1:nLocalSteps
    movBody = buildBodyStruct(Fmov, VmovP);
    movBody.bbox = bboxFromV(VmovP);

    % ---- Run contact area ----
    res = computeContactArea_STS(cubeBodyPlaced, movBody, tol, opts);

    A = NaN;
    mask = [];

    if isstruct(res)
        if isfield(res,'contactArea'), A = res.contactArea; end

        % Try common mask field names (your code may differ)
        if isfield(res,'contactMask'), mask = res.contactMask; end
        if isempty(mask) && isfield(res,'maskSlave'), mask = res.maskSlave; end
        if isempty(mask) && isfield(res,'mask'), mask = res.mask; end
    else
        A = res;
    end

    % ---- Fallback mask if none returned ----
    if isempty(mask)
        warning('computeContactArea_STS did not return a contact mask. Using distance-based fallback mask for visualisation/refinement.');
        mask = fallbackMaskByDistance(cubeBodyPlaced, movBody, tol, fallbackDistFactor);
    end

    fprintf('Step %d/%d | faces=%d | A=%.4f mm^2 | contactTris=%d\n', ...
        s, nLocalSteps, size(Fmov,1), A, nnz(mask));

    snap(s).F = Fmov;
    snap(s).V = VmovP;
    snap(s).mask = mask;
    snap(s).A = A;

    % Stop if last step (don't refine further)
    if s == nLocalSteps
        break
    end

    % If no contact triangles, refinement can't proceed
    if ~any(mask)
        warning('No contact triangles detected at step %d. Stopping refinement.', s);
        break
    end

    % ---- Refine ONLY those triangles ----
    [Fmov, VmovP] = refineTrianglesMidpoint_cached(Fmov, VmovP, mask);
end

%% ===== PLOTS: coarse, step-before-final, final =====
n = numel(snap);
if n == 0
    error('No steps recorded. Something failed before first contact run.');
end

plotStep(snap, 1,            cubeBodyPlaced, tol, 'COARSE (step 1)');
if n >= 2
    plotStep(snap, max(1,n-1), cubeBodyPlaced, tol, 'STEP BEFORE FINAL');
end
plotStep(snap, n,            cubeBodyPlaced, tol, 'FINAL');

fprintf('\nDone.\n');

%% ========================= FUNCTIONS =========================

function [F,V] = loadMeshAny(path)
    if exist('loadStlMesh','file')
        [F,V] = loadStlMesh(path);
    elseif exist('loadAnyStl','file')
        [F,V] = loadAnyStl(path);
    else
        error('Neither loadStlMesh nor loadAnyStl found on path.');
    end
end

function [Fsurf, Vsurf] = surfaceFromStl_generateMesh(stlPath, Hmax)
    model = createpde();
    importGeometry(model, stlPath);

    try
        msh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'linear');
    catch
        msh = generateMesh(model, 'Hmax', Hmax);
    end

    V = msh.Nodes.';           % Nx3
    E = msh.Elements;          % nNodesPerElem x nElems
    nper = size(E,1);

    if nper == 4
        T = E.';               % linear tets
    elseif nper > 4
        T = E(1:4,:).';        % corner nodes
    elseif nper == 3
        Fsurf = E.'; Vsurf = V;
        return
    else
        error('Unexpected msh.Elements size: %dx%d', size(E,1), size(E,2));
    end

    TR = triangulation(T, V);
    [Fsurf, Vsurf] = freeBoundary(TR);
end

function [VcubeP, VsphP] = placeSphereOnCubeTop(Vcube, Vsph, penetration)
% Align XY centroids, then sit sphere on cube top with penetration.
    VcubeP = Vcube;

    cC = mean(Vcube(:,1:2),1);
    cS = mean(Vsph(:,1:2),1);

    VsphP = Vsph;
    VsphP(:,1:2) = VsphP(:,1:2) + (cC - cS);

    zTopCube   = max(VcubeP(:,3));
    zBottomSph = min(VsphP(:,3));

    VsphP(:,3) = VsphP(:,3) + ((zTopCube - penetration) - zBottomSph);
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
end

function mask = fallbackMaskByDistance(masterBody, movBody, tol, distFactor)
% Approximate "contact triangles" by centroid distance to nearest master centroid.
    if ~isfield(masterBody,'kdtree') || isempty(masterBody.kdtree)
        masterBody.kdtree = KDTreeSearcher(masterBody.centroids);
    end
    C = movBody.centroids;
    [~, d] = knnsearch(masterBody.kdtree, C, 'K', 1);
    mask = d <= distFactor * tol;
end

function [F2, V2] = refineTrianglesMidpoint_cached(F, V, refineMask)
% Refine ONLY triangles in refineMask using midpoint subdivision.
% Uses an edge-midpoint cache so shared edges re-use the same midpoint vertex.

    idxRef  = find(refineMask);
    idxKeep = find(~refineMask);

    Fkeep = F(idxKeep,:);
    Fr = F(idxRef,:);
    V2 = V;

    % Map edge key -> midpoint vertex index
    edgeMap = containers.Map('KeyType','char','ValueType','int32');

    Fnew = zeros(4*size(Fr,1), 3);
    outRow = 0;

    for i = 1:size(Fr,1)
        a = Fr(i,1); b = Fr(i,2); c = Fr(i,3);

        iab = getMid(a,b);
        ibc = getMid(b,c);
        ica = getMid(c,a);

        % 4 child triangles
        outRow = outRow + 1; Fnew(outRow,:) = [a,   iab, ica];
        outRow = outRow + 1; Fnew(outRow,:) = [iab, b,   ibc];
        outRow = outRow + 1; Fnew(outRow,:) = [ica, ibc, c];
        outRow = outRow + 1; Fnew(outRow,:) = [iab, ibc, ica];
    end

    F2 = [Fkeep; Fnew(1:outRow,:)];

    % remove unused vertices (keeps arrays tidy)
    [V2, F2] = removeUnusedVertices(V2, F2);

    function midIdx = getMid(i1,i2)
        if i1 < i2
            key = sprintf('%d_%d', i1, i2);
        else
            key = sprintf('%d_%d', i2, i1);
        end

        if isKey(edgeMap, key)
            midIdx = edgeMap(key);
            return
        end

        vm = 0.5*(V2(i1,:) + V2(i2,:));
        midIdx = size(V2,1) + 1;
        V2(midIdx,:) = vm;
        edgeMap(key) = midIdx;
    end
end

function [V2, F2] = removeUnusedVertices(V, F)
    used = false(size(V,1),1);
    used(F(:)) = true;
    map = zeros(size(V,1),1);
    map(used) = 1:nnz(used);
    V2 = V(used,:);
    F2 = map(F);
end

function plotStep(snap, idx, masterBody, tol, labelStr)
    Fm = snap(idx).F;
    Vm = snap(idx).V;
    mask = snap(idx).mask;
    A = snap(idx).A;

    % Decimate for display only (prevents rendering issues)
    [Fmd, Vmd] = decimateForPlot(Fm, Vm, 50000);
    mask_d = mapMaskToDecimated(Fm, Fmd, mask); % best-effort mapping

    fig = figure('Color','w');
    set(fig,'Name',sprintf('Sphere–Cube | %s', labelStr));
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | tol=%.3f mm | A=%.4f mm^2 | step=%d', labelStr, tol, A, idx), 'Interpreter','none');

    % Master cube (faint)
    [Fcd, Vcd] = decimateForPlot(masterBody.F, masterBody.V, 60000);
    patch('Faces', Fcd, 'Vertices', Vcd, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.12);

    % Moving sphere (faint green)
    patch('Faces', Fmd, 'Vertices', Vmd, ...
        'FaceColor', [0 0.7 0], 'EdgeColor', 'none', 'FaceAlpha', 0.12);

    % Highlight "contact triangles" (magenta)
    if any(mask_d)
        idxTri = find(mask_d);
        patch('Faces', Fmd(idxTri,:), 'Vertices', Vmd, ...
            'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end

    view(3); camlight; lighting gouraud
    drawnow
end

function [F2,V2] = decimateForPlot(F,V,targetFaces)
    if size(F,1) <= targetFaces
        F2 = F; V2 = V; return
    end
    p = targetFaces / size(F,1);
    [F2,V2] = reducepatch(F,V,p);
end

function mask_d = mapMaskToDecimated(Ffull, Fdec, mask_full)
% Best-effort mapping of triangle mask after reducepatch.
% reducepatch changes faces, so exact mapping isn't guaranteed.
% For display this is OK: we re-classify triangles in decimated mesh by
% checking if their vertices existed in any "contact triangle" vertex set.

    if isempty(mask_full) || ~any(mask_full)
        mask_d = false(size(Fdec,1),1);
        return
    end

    contactVerts = unique(Ffull(mask_full,:));
    mask_d = ismember(Fdec(:,1), contactVerts) & ismember(Fdec(:,2), contactVerts) & ismember(Fdec(:,3), contactVerts);
end
