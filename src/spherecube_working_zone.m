% sphereCube_localRefine_penetrationArea_autoStop.m
% ------------------------------------------------------------
% Sphere–Cube local refinement using PENETRATION CONTACT AREA:
% - Contact area = area of the sphere surface that lies "inside" the cube
% - NEW method: per-triangle HALF-SPACE CLIPPING (partial triangle areas)
% - Refinement: refine ONLY the triangles that contributed contact area
% - Auto-stop: dA% < targetPct for needConsecutive steps
%
% Outputs:
%  1) Per-step mesh/contact figure (optional)
%  2) Coarse contact
%  3) Step-before-final contact
%  4) Final contact
%  5) Final zoned refinement plot (faces coloured by refinement level)
%
% REQUIREMENTS:
%  - buildBodyStruct
%  - loadStlMesh OR loadAnyStl
% ------------------------------------------------------------

clear; clc; close all;

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')));
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Sphere–Cube LOCAL REFINEMENT (PENETRATION AREA, auto-stop) ===\n');
fprintf('Project root: %s\n',projectRoot);
fprintf('Model dir:    %s\n\n',modelDir);

%% ===== FILES =====
cubeStl   = fullfile(modelDir,'50mm cube 4 iterations.stl');      % MASTER
sphereStl = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');   % MOVING

assert(exist(cubeStl,'file')==2,'Cube STL missing: %s', cubeStl);
assert(exist(sphereStl,'file')==2,'Sphere STL missing: %s', sphereStl);

%% ===== PENETRATION SETUP =====
% We still need to place the sphere so it slightly penetrates the cube.
penetration = 5;   % mm (set whatever you want to test)
fprintf('Using penetration = %.3f mm\n', penetration);

%% ===== LOCAL REFINEMENT CONTROL =====
HmaxCoarse  = 6.0;          % mm initial sphere mesh (generateMesh surface)

% ---- AUTO STOP SETTINGS ----
maxSteps        = 25;       % safety cap
targetPct       = 5.0;      % stop when % change in area < 5%
needConsecutive = 4;        % require this many consecutive steps under target
minValidArea    = 1e-12;    % ignore 0/near-0 areas in convergence check

%% ===== (OPTIONAL) PER-STEP DISPLAY SETTINGS =====
showEveryStep   = true;     % show per-step figures
pauseEachFigure = false;    % pause each step
closePrevFigure = true;     % keep only one per-step window open
stepFigHandle   = [];

%% ===== LOAD MASTER CUBE =====
if ~exist('buildBodyStruct','file'), error('buildBodyStruct not found on path.'); end

[Fcube,Vcube] = loadMeshAny(cubeStl);
cubeBody = buildBodyStruct(Fcube,Vcube);

% Precompute normals + KD-tree on cube triangle centroids (for "nearest plane")
cubeBody.F = Fcube;
cubeBody.V = Vcube;
cubeBody.faceNormals = faceNormalsFromFV(Fcube,Vcube);
cubeBody.kdtreeFaces = KDTreeSearcher(cubeBody.centroids);

%% ===== INITIAL COARSE SPHERE SURFACE =====
[Fmov,Vmov] = surfaceFromStl_generateMesh(sphereStl,HmaxCoarse);

% refinement level per face (0 = coarse)
faceLvl = zeros(size(Fmov,1),1,'uint16');

% Place sphere on cube top with specified penetration
[VcubeP,VmovP] = placeSphereOnCubeTop(Vcube,Vmov,penetration);
cubeBody.V = VcubeP;  % update master vertices after placement

%% ===== ITERATIVE: PENETRATION AREA → LOCAL REFINE (AUTO STOP) =====
snap = struct('F',{},'V',{},'mask',{},'A',{},'faceLvl',{},'dPct',{});

belowCount = 0;
Aprev = NaN;

for s = 1:maxSteps
    % --- penetration contact area on current moving mesh ---
    [A, mask] = computePenetrationArea_clipHalfSpace(cubeBody, Fmov, VmovP);

    % --- % change vs previous valid area ---
    dPct = NaN;
    if ~isnan(Aprev) && ~isnan(A) && Aprev > minValidArea
        dPct = 100 * abs(A - Aprev) / Aprev;
    end

    fprintf('Step %d/%d | faces=%d | contactFaces=%d | Apen=%.6f mm^2 | dA=%.2f%%\n', ...
        s, maxSteps, size(Fmov,1), nnz(mask), A, dPct);

    % --- store snapshot ---
    snap(s).F = Fmov;
    snap(s).V = VmovP;
    snap(s).mask = mask;
    snap(s).A = A;
    snap(s).faceLvl = faceLvl;
    snap(s).dPct = dPct;

    % --- per-step plot ---
    if showEveryStep
        if closePrevFigure && ~isempty(stepFigHandle) && isvalid(stepFigHandle)
            close(stepFigHandle);
        end
        ttl = sprintf('STEP %d | faces=%d | contactFaces=%d | Apen=%.6f | dA=%.2f%%', ...
            s, size(Fmov,1), nnz(mask), A, dPct);
        stepFigHandle = plotStepFigure(Fmov, VmovP, mask, cubeBody, ttl);
        drawnow;
        if pauseEachFigure, disp('Press any key to continue...'); pause; end
    end

    % --- stop if no penetration faces detected ---
    if ~any(mask)
        warning('No penetration triangles detected at step %d. Stopping.', s);
        break
    end

    % --- stopping rule: below targetPct for needConsecutive steps ---
    if ~isnan(dPct) && dPct < targetPct
        belowCount = belowCount + 1;
    else
        belowCount = 0;
    end

    if ~isnan(A) && A > minValidArea
        Aprev = A;
    end

    if belowCount >= needConsecutive
        fprintf('STOP: dA < %.2f%% for %d consecutive steps.\n', targetPct, needConsecutive);
        break
    end

    % --- refine ONLY faces that contributed penetration area ---
    [Fmov, VmovP, faceLvl] = refineTrianglesMidpoint_cached(Fmov, VmovP, mask, faceLvl);
end

% Trim snap if loop ended early
snap = snap(1:min(numel(snap), s));

%% ===== SUMMARY VISUALS =====
plotStepSummary(snap, 1, cubeBody, 'COARSE (PENETRATION AREA)');
nSnap = numel(snap);
if nSnap >= 2
    plotStepSummary(snap, max(1,nSnap-1), cubeBody, 'STEP BEFORE FINAL');
end
plotStepSummary(snap, nSnap, cubeBody, 'FINAL (PENETRATION AREA)');
plotFinalZoning(snap(nSnap), cubeBody);

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
        Fsurf = E.'; Vsurf = V; return
    end
    T = E(1:4,:).';
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

function N = faceNormalsFromFV(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    N = cross(v2 - v1, v3 - v1, 2);
    nrm = sqrt(sum(N.^2,2));
    nrm(nrm==0) = 1;
    N = N ./ nrm;
end

% -------------------------------------------------------------------------
% NEW METHOD: Penetration contact area via per-triangle half-space clipping
% -------------------------------------------------------------------------
function [Apen, maskPen] = computePenetrationArea_clipHalfSpace(master, Fslave, Vslave)
% For each slave triangle:
%  1) find nearest master face (by centroid KD-tree)
%  2) use that face plane (point + normal) as local separating plane
%  3) compute signed distances of the 3 vertices
%  4) clip triangle by d<=0 and add clipped area
%
% Returns:
%   Apen     penetration contact area (mm^2)
%   maskPen  faces that contributed non-zero clipped area

    % Slave triangles
    P1 = Vslave(Fslave(:,1),:);
    P2 = Vslave(Fslave(:,2),:);
    P3 = Vslave(Fslave(:,3),:);

    C  = (P1 + P2 + P3) / 3;  % slave face centroids

    % Nearest master face for each slave face centroid
    idx = knnsearch(master.kdtreeFaces, C, 'K', 1);

    % master plane: point on plane + normal
    Fm = master.F(idx,:);
    V1m = master.V(Fm(:,1),:);
    Nm  = master.faceNormals(idx,:);  % assumed outward

    % Signed distances for each slave vertex to that local plane
    d1 = dot(P1 - V1m, Nm, 2);
    d2 = dot(P2 - V1m, Nm, 2);
    d3 = dot(P3 - V1m, Nm, 2);

    Aclip = zeros(size(Fslave,1),1);
    for i = 1:size(Fslave,1)
        Aclip(i) = clippedAreaTriangleHalfSpace(P1(i,:),P2(i,:),P3(i,:), d1(i),d2(i),d3(i));
    end

    Apen = sum(Aclip);
    maskPen = Aclip > 0;
end

function A = clippedAreaTriangleHalfSpace(p1,p2,p3, d1,d2,d3)
% Clip triangle by half-space d<=0, where d are signed distances at vertices
% Linear interpolation along edges.
% Result polygon has 0..4 vertices; compute its area in 3D.

    pts = [p1; p2; p3];
    ds  = [d1; d2; d3];

    inside = ds <= 0;

    if all(~inside)
        A = 0; return
    elseif all(inside)
        A = triArea3D(p1,p2,p3); return
    end

    poly = [];  % clipped polygon vertices in 3D
    for k = 1:3
        k2 = mod(k,3) + 1;
        Pk  = pts(k,:);  dk  = ds(k);
        Pk2 = pts(k2,:); dk2 = ds(k2);

        if dk <= 0
            poly = [poly; Pk];
        end

        % Edge crosses boundary?
        if (dk <= 0 && dk2 > 0) || (dk > 0 && dk2 <= 0)
            t = dk / (dk - dk2);  % where d(t)=0
            Pi = Pk + t*(Pk2 - Pk);
            poly = [poly; Pi];
        end
    end

    % Polygon area: triangulate fan around poly(1)
    if size(poly,1) < 3
        A = 0; return
    end
    A = 0;
    p0 = poly(1,:);
    for j = 2:(size(poly,1)-1)
        A = A + triArea3D(p0, poly(j,:), poly(j+1,:));
    end
end

function A = triArea3D(a,b,c)
    A = 0.5 * norm(cross(b-a, c-a));
end

% -------------------------------------------------------------------------
% Refinement (unchanged)
% -------------------------------------------------------------------------
function [F2,V2,lvl2] = refineTrianglesMidpoint_cached(F,V,mask,lvl)
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
    lvl2 = lvl;
end

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
function h = plotStepFigure(F,V,mask,master,ttl)
    h = figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(ttl,'Interpreter','none');

    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    patch('Faces', F, 'Vertices', V, ...
        'FaceColor', [0 .7 0], 'FaceAlpha', 0.06, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.22, 'LineWidth', 0.25);

    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], 'FaceAlpha', 0.92, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 0.25, 'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end

function plotStepSummary(snap,idx,master,label)
    F = snap(idx).F; V = snap(idx).V; mask = snap(idx).mask; A = snap(idx).A;

    figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | Apen=%.6f mm^2', label, A), 'Interpreter','none');

    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    patch('Faces', F, 'Vertices', V, ...
        'FaceColor', [0 .7 0], 'FaceAlpha', 0.08, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, 'LineWidth', 0.25);

    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], 'FaceAlpha', 0.90, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 0.25, 'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end

function plotFinalZoning(snap,master)
    F = snap.F; V = snap.V; lvl = double(snap.faceLvl); mask = snap.mask; A = snap.A;

    figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('FINAL ZONED LOCAL REFINEMENT | Apen=%.6f mm^2', A), 'Interpreter','none');

    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    patch('Faces', F, 'Vertices', V, ...
        'FaceVertexCData', lvl, 'FaceColor', 'flat', ...
        'FaceAlpha', 0.35, 'EdgeColor', [0 0 0], ...
        'EdgeAlpha', 0.22, 'LineWidth', 0.25);

    cb = colorbar;
    cb.Label.String = 'Refinement level (0=coarse, higher=finer)';

    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], 'FaceAlpha', 0.95, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 0.25, 'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end
