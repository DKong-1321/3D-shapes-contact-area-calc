% Spheregroove.m
% -------------------------------------------------------------------------
% Mesh convergence for ONE PAIR (sphere in hemispherical hole) using:
%   - generateMesh ONCE to get an initial coarse MOVING sphere surface
%   - compute PENETRATION "contact area" via half-space clipping per-triangle
%   - refine ONLY triangles that contribute to penetration area
%   - stop when dA% < convTolPct for convConsec consecutive levels (after minLevels)
%
% Placement:
%   - rotate fixed groove about X
%   - align sphere in XY to cavity region
%   - place sphere bottom (min Z) at cavity-bottom (local min Z), not global structure min Z
%   - optionally push in with zEps < 0 to create penetration patch
%
% Plots:
%   - debug plot EVERY level (wireframes)
%   - COARSE / STEP BEFORE FINAL / FINAL contact-highlight plot
%   - convergence plots using mean edge length of CONTACT triangles only
%   - FINAL ZONING plot (refinement levels + last contact patch)
%
% REQUIREMENTS:
%   - buildBodyStruct (must exist in your /src folder)
%   - loadStlMesh OR loadAnyStl (must exist in your /src folder)
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% ======================= PATHS =======================
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Mesh convergence | sphere in hemi-hole ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

%% ======================= FILES =======================
sphereStlFine = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');      % MOVING
hemiHoleStl   = fullfile(modelDir,['hemisphere groove 50mm dia' ....stl']);      % FIXED

assert(exist(sphereStlFine,'file')==2, 'Sphere STL not found: %s', sphereStlFine);
assert(exist(hemiHoleStl,'file')==2,   'Hemi-hole STL not found: %s', hemiHoleStl);

fprintf('STLs:\n  Sphere:    %s\n  Hemi-hole: %s\n\n', sphereStlFine, hemiHoleStl);

%% ======================= SETTINGS =======================
tol == 0.20;          % (label only; penetration method doesn't use tol band)
zEps == -0.2;         % mm. 0 means just-touch bottom; negative pushes INTO cavity for a patch.

opts = struct();
opts.tol = tol;
opts.debugPlot = true;
opts.debugPlotEveryLevel = true;
opts.debugPlotPause = false;
opts.showRoiBox = false;

HmaxStart  = 6.0;    % initial sphere mesh coarseness
maxLevels  = 25;     % safety cap
convTolPct = 5;      % <5% change
convConsec = 3;      % 3 in a row
minLevels  = 3;

fprintf('Settings:\n');
fprintf('  HmaxStart=%.3f mm | maxLevels=%d | conv=%.2f%% for %d in a row | minLevels=%d\n', ...
    HmaxStart, maxLevels, convTolPct, convConsec, minLevels);
fprintf('  zEps = %.3f mm (negative pushes sphere down into cavity)\n\n', zEps);

%% ======================= LOAD FIXED BODY =======================
[Fh,Vh] = loadMeshAny(hemiHoleStl);

% Rotate groove about X so cavity faces up (flip sign if needed)
Vh = rotateAboutCentroid(Vh, rotxd(90));

hemiMaster0 = buildBodyStruct(Fh,Vh);

%% ======================= RUN PAIR =======================
fprintf('--- Pair: sphere-in-hemiHole ---\n');

Result = runPair_refineContactTris_only( ...
    'sphere-in-hemiHole', ...
    hemiMaster0, ...
    sphereStlFine, ...
    @(Vfix,Vmov) placeSphereInHemiHole_localBottomAlign(Vfix,Vmov,zEps), ...
    HmaxStart, ...
    maxLevels, ...
    convTolPct, ...
    convConsec, ...
    minLevels, ...
    tol, ...
    opts, ...
    2100 ...
);

%% ======================= CONVERGENCE PLOTS =======================
figure('Color','w'); hold on; grid on
xlabel('log_{10}(mean edge length of CONTACT triangles)  [mm]');
ylabel('Penetration contact area (mm^2)');
title(sprintf('Mesh convergence | sphere-in-hemiHole | tol(label)=%.3f mm', tol));
plot(log10(Result.meanEdgeLen_contact_mm), Result.A, '-o');
set(gca,'XDir','reverse');

figure('Color','w'); hold on; grid on
xlabel('log_{10}(mean edge length of CONTACT triangles)  [mm]');
ylabel('% change vs previous level');
title('Convergence rate (% change)');
plot(log10(Result.meanEdgeLen_contact_mm), Result.dA_pct, '-o');
set(gca,'XDir','reverse');

%% ======================= FINAL ZONING =======================
plotFinalZoningPair(3001, Result, tol);

%% ======================= SAVE =======================
outPath = fullfile(projectRoot,'results_meshConvergence_sphereInHemiHole_contactEdgeOnly.mat');
save(outPath,'Result','tol','HmaxStart','convTolPct','convConsec','minLevels','maxLevels','zEps','sphereStlFine','hemiHoleStl');
fprintf('\nSaved: %s\n', outPath);
fprintf('=== DONE ===\n');

%% =================================================================
%% ============================ FUNCTIONS ===========================
%% =================================================================

function [F,V] = loadMeshAny(path)
    if exist('loadStlMesh','file')
        [F,V] = loadStlMesh(path);
    elseif exist('loadAnyStl','file')
        [F,V] = loadAnyStl(path);
    else
        error('Neither loadStlMesh nor loadAnyStl found on path.');
    end
end

function out = runPair_refineContactTris_only(name, masterBody0, movStlPath, placer, ...
    HmaxStart, maxLevels, convTolPct, convConsec, minLevels, tol, opts, figBase)

    Vfix0 = masterBody0.V;

    A   = nan(1,maxLevels);
    Nf  = nan(1,maxLevels);
    dA  = nan(1,maxLevels);

    EmC  = nan(1,maxLevels);
    E50C = nan(1,maxLevels);
    EmAll = nan(1,maxLevels);

    % initial coarse moving surface (UNPLACED)
    [Fmov, Vmov] = surfaceFromStl_generateMesh(movStlPath, HmaxStart);
    faceLvl = zeros(size(Fmov,1),1,'uint16');

    % PLACE ONCE, FREEZE TRANSFORM
    [VfixP0, VmovP0] = placer(Vfix0, Vmov);

    % Freeze moving translation using bbox-centre (reduces drift)
    tMove = bboxCenter(VmovP0) - bboxCenter(Vmov);
    VfixPlacedConst = VfixP0;

    V_masterPlaced = VfixPlacedConst;
    V_movPlaced    = VmovP0;

    snap = struct('F',{},'V',{},'mask',{},'A',{},'faceLvl',{},'dPct',{});
    belowCount = 0;
    Aprev = NaN;
    mask_last = [];

    for k = 1:maxLevels
        % Apply frozen placement
        VfixP = VfixPlacedConst;
        VmovP = Vmov + tMove;

        masterBody = buildBodyStruct(masterBody0.F, VfixP);
        masterBody.bbox = bboxFromV(VfixP);
        masterBody.faceNormals = faceNormalsFromFV(masterBody.F, masterBody.V);
        masterBody.kdtreeFaces = KDTreeSearcher(masterBody.centroids);

        movBody = buildBodyStruct(Fmov, VmovP);
        movBody.bbox = bboxFromV(VmovP);

        V_masterPlaced = VfixP;
        V_movPlaced    = VmovP;

        Nf(k) = size(Fmov,1);

        if isfield(opts,'debugPlot') && opts.debugPlot ...
                && isfield(opts,'debugPlotEveryLevel') && opts.debugPlotEveryLevel
            debugPlotPair(figBase + k, name, k, masterBody, movBody, tol, opts);
        end

        % Penetration area + mask
        res = computeContactArea_STS_penetration(masterBody, movBody, tol, opts);
        mask = res.contactMask;
        A(k) = res.contactArea;

        if isempty(mask) || ~islogical(mask) || numel(mask) ~= size(Fmov,1)
            mask = false(size(Fmov,1),1);
            warning('%s: contact mask invalid at level %d.', name, k);
        end
        mask_last = mask;

        % Edge stats on CONTACT faces only
        if any(mask)
            [EmC(k), E50C(k)] = meshEdgeStats_faces(Fmov(mask,:), VmovP);
        else
            EmC(k)  = NaN;
            E50C(k) = NaN;
        end
        EmAll(k) = meshEdgeStats_all(Fmov, VmovP);

        % Convergence metric
        dPct = NaN;
        if ~isnan(Aprev) && ~isnan(A(k)) && Aprev > 1e-12
            dPct = 100 * abs(A(k) - Aprev) / Aprev;
        end
        dA(k) = dPct;

        fprintf('  Level %d | faces=%d | meanEdge(CONTACT)=%.3f mm | Apen=%.6f | dA%%=%.3f\n', ...
            k, Nf(k), EmC(k), A(k), dPct);

        snap(k).F = Fmov;
        snap(k).V = VmovP;
        snap(k).mask = mask;
        snap(k).A = A(k);
        snap(k).faceLvl = faceLvl;
        snap(k).dPct = dPct;

        if ~isnan(A(k)) && A(k) > 1e-12
            Aprev = A(k);
        end

        if ~any(mask)
            warning('%s: No penetration triangles detected at level %d. Stopping.', name, k);
            A = A(1:k); Nf = Nf(1:k); dA = dA(1:k);
            EmC = EmC(1:k); E50C = E50C(1:k); EmAll = EmAll(1:k);
            snap = snap(1:k);
            break
        end

        % Convergence counter
        if k >= minLevels && ~isnan(dPct) && dPct < convTolPct
            belowCount = belowCount + 1;
        else
            belowCount = 0;
        end

        if k >= minLevels && belowCount >= convConsec
            fprintf('  -> Converged: dA%% < %.2f for %d consecutive levels. Stopping.\n', convTolPct, convConsec);
            A = A(1:k); Nf = Nf(1:k); dA = dA(1:k);
            EmC = EmC(1:k); E50C = E50C(1:k); EmAll = EmAll(1:k);
            snap = snap(1:k);
            break
        end

        % Refine ONLY contributing faces (UNPLACED mesh)
        [Fmov, Vmov, faceLvl] = refineTrianglesMidpoint_cached(Fmov, Vmov, mask, faceLvl);
    end

    nSnap = numel(snap);

    % show COARSE / STEP BEFORE FINAL / FINAL
    plotStepSphereCubeStyle(snap, 1, masterBody0.F, V_masterPlaced, tol, sprintf('%s | COARSE', name));
    if nSnap >= 2
        plotStepSphereCubeStyle(snap, max(1,nSnap-1), masterBody0.F, V_masterPlaced, tol, sprintf('%s | STEP BEFORE FINAL', name));
    end
    plotStepSphereCubeStyle(snap, nSnap, masterBody0.F, V_masterPlaced, tol, sprintf('%s | FINAL CONTACT', name));

    out = struct('name',name, ...
        'A',A(1:nSnap), ...
        'Nfaces_moving',Nf(1:nSnap), ...
        'dA_pct',dA(1:nSnap), ...
        'meanEdgeLen_contact_mm',EmC(1:nSnap), ...
        'medianEdgeLen_contact_mm',E50C(1:nSnap), ...
        'meanEdgeLen_all_mm',EmAll(1:nSnap), ...
        'F_final',Fmov, ...
        'V_finalPlaced',V_movPlaced, ...
        'faceLvl_final',faceLvl, ...
        'F_master',masterBody0.F, ...
        'V_masterPlaced',V_masterPlaced, ...
        'mask_last',mask_last);
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
        T = E(1:4,:).';        % corner nodes only
    elseif nper == 3
        Fsurf = E.'; Vsurf = V; return
    else
        error('Unexpected msh.Elements size: %dx%d', size(E,1), size(E,2));
    end

    TR = triangulation(T, V);
    [Fsurf, Vsurf] = freeBoundary(TR);
end

function [VfixP, VmovP] = placeSphereInHemiHole_localBottomAlign(Vfix, Vmov, zEps)
% Align sphere to cavity region; put sphere bottom at local cavity bottom.
    VfixP = Vfix;

    % sphere radius estimate
    bbS = bboxFromV(Vmov);
    Rs  = 0.5 * (bbS(1,2) - bbS(1,1));

    % rough cavity XY centre from fixed bbox
    c0 = bboxCenter(VfixP);

    % choose candidate cavity region: lower part + near bbox centre
    zCut = prctile(VfixP(:,3), 60);
    idxLow = VfixP(:,3) <= zCut;

    dxy0 = hypot(VfixP(:,1)-c0(1), VfixP(:,2)-c0(2));
    idxNear = dxy0 <= (1.5*Rs);

    idx = idxLow & idxNear;
    if nnz(idx) > 50
        cXY = mean(VfixP(idx,1:2), 1);
    else
        cXY = c0(1:2);
    end

    % Align sphere XY
    cMov = bboxCenter(Vmov);
    VmovP = Vmov;
    VmovP(:,1:2) = VmovP(:,1:2) + (cXY - cMov(1:2));

    % Local cavity bottom Z from vertices near cavity footprint
    dxy = hypot(VfixP(:,1)-cXY(1), VfixP(:,2)-cXY(2));
    idxLocal = dxy <= (1.05*Rs);

    if nnz(idxLocal) > 50
        zCavBottom = min(VfixP(idxLocal,3));
    else
        zCavBottom = min(VfixP(:,3));
        warning('Local cavity selection small; falling back to global min Z.');
    end

    % Bottom-align sphere
    zSphereBottom = min(VmovP(:,3));
    VmovP(:,3) = VmovP(:,3) + ((zCavBottom + zEps) - zSphereBottom);
end

function meanE = meshEdgeStats_all(F, V)
    if isempty(F) || isempty(V), meanE = NaN; return; end
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2,2));
    meanE = mean(L);
end

function [meanE, medE] = meshEdgeStats_faces(Fsub, V)
    if isempty(Fsub) || isempty(V), meanE = NaN; medE = NaN; return; end
    E = [Fsub(:,[1 2]); Fsub(:,[2 3]); Fsub(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2,2));
    meanE = mean(L);
    medE  = median(L);
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
end

function c = bboxCenter(V)
    bb = bboxFromV(V);
    c  = [mean(bb(1,:)), mean(bb(2,:)), mean(bb(3,:))];
end

function debugPlotPair(figId, pairName, levelIdx, fixedBody, movBody, tol, opts)
    figure(figId); clf
    set(gcf,'Name',sprintf('%s | level %d', pairName, levelIdx));
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | level %d | tol(label)=%.3f', pairName, levelIdx, tol), 'Interpreter','none');

    patch('Faces', fixedBody.F, 'Vertices', fixedBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.35);

    patch('Faces', movBody.F, 'Vertices', movBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0.7 0], 'EdgeAlpha', 0.35);

    view(3);
end

function R = rotxd(deg)
    a = deg2rad(deg);
    R = [1 0 0;
         0 cos(a) -sin(a);
         0 sin(a)  cos(a)];
end

function V2 = rotateAboutCentroid(V, R)
    c  = mean(V,1);
    V2 = (R * (V - c)')' + c;
end

function [F2,V2,lvl2] = refineTrianglesMidpoint_cached(F,V,mask,lvl)
    if nargin < 4 || isempty(lvl)
        lvl = zeros(size(F,1),1,'uint16');
    end

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
            m = edgeMap(k); return
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

function plotStepSphereCubeStyle(snap,idx,Fmaster,VmasterPlaced,tol,label)
    F = snap(idx).F;
    V = snap(idx).V;
    mask = snap(idx).mask;
    A = snap(idx).A;

    figure('Color','w'); hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | tol(label)=%.3f mm | A=%.6f mm^2', label, tol, A), 'Interpreter','none');

    patch('Faces', Fmaster, 'Vertices', VmasterPlaced, ...
        'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    patch('Faces', F, 'Vertices', V, ...
        'FaceColor', [0 .7 0], 'FaceAlpha', 0.08, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, ...
        'LineWidth', 0.25);

    if any(mask)
        patch('Faces', F(mask,:), 'Vertices', V, ...
            'FaceColor', [1 0 1], 'FaceAlpha', 0.90, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 0.25, ...
            'LineWidth', 0.30);
    end

    view(3); camlight; lighting gouraud
end

function plotFinalZoningPair(figId, R, tol)
    figure(figId); clf
    set(gcf,'Color','w','Name',['ZONED: ',R.name]);
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | FINAL zoned refinement | tol(label)=%.3f mm', R.name, tol), 'Interpreter','none');

    patch('Faces',R.F_master,'Vertices',R.V_masterPlaced, ...
        'FaceColor','none','EdgeColor',[0 0 0],'EdgeAlpha',0.18,'LineWidth',0.25);

    F = R.F_final;
    V = R.V_finalPlaced;
    lvl = double(R.faceLvl_final);

    patch('Faces',F,'Vertices',V, ...
        'FaceVertexCData',lvl, ...
        'FaceColor','flat', ...
        'EdgeColor',[0 0 0], ...
        'EdgeAlpha',0.25, ...
        'FaceAlpha',0.35, ...
        'LineWidth',0.25);

    cb = colorbar;
    cb.Label.String = 'Refinement level (0=coarse, higher=finer)';

    if isfield(R,'mask_last') && ~isempty(R.mask_last) && any(R.mask_last)
        patch('Faces',F(R.mask_last,:), 'Vertices', V, ...
            'FaceColor',[1 0 1], 'EdgeColor','none', 'FaceAlpha',0.90);
    end

    view(3); camlight; lighting gouraud
end

%% ===================== CONTACT FUNCTION (CLIPPED PENETRATION AREA) =====================

function results = computeContactArea_STS_penetration(master, slave, varargin)
    tol  = [];
    opts = struct();

    if numel(varargin) == 1
        if isstruct(varargin{1})
            opts = varargin{1};
        elseif isnumeric(varargin{1}) && isscalar(varargin{1})
            tol = varargin{1};
        else
            error('Invalid third argument. Expected opts struct or scalar tol.');
        end
    elseif numel(varargin) == 2
        if ~(isnumeric(varargin{1}) && isscalar(varargin{1}))
            error('If 4 inputs are used, the 3rd must be scalar tol.');
        end
        tol = varargin{1};
        if isstruct(varargin{2})
            opts = varargin{2};
        else
            error('If 4 inputs are used, the 4th must be opts struct.');
        end
    elseif numel(varargin) > 2
        error('Too many input arguments.');
    end

    if ~isempty(tol)
        opts.tol = tol; %#ok<NASGU>
    end

    master = ensureBodyFields_fast(master);
    slave  = ensureBodyFields_fast(slave);

    if ~isfield(master,'faceNormals') || isempty(master.faceNormals)
        master.faceNormals = faceNormalsFromFV(master.F, master.V);
    end
    if ~isfield(master,'kdtreeFaces') || isempty(master.kdtreeFaces)
        master.kdtreeFaces = KDTreeSearcher(master.centroids);
    end

    [Apen, mask] = computePenetrationArea_clipHalfSpace(master, slave.F, slave.V);

    results = struct();
    results.contactArea = Apen;
    results.contactMask = mask;
end

function [Apen, maskPen] = computePenetrationArea_clipHalfSpace(master, Fslave, Vslave)
    P1 = Vslave(Fslave(:,1),:);
    P2 = Vslave(Fslave(:,2),:);
    P3 = Vslave(Fslave(:,3),:);
    C  = (P1 + P2 + P3) / 3;

    idx = knnsearch(master.kdtreeFaces, C, 'K', 1);

    Fm  = master.F(idx,:);
    V1m = master.V(Fm(:,1),:);
    Nm  = master.faceNormals(idx,:);

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

function N = faceNormalsFromFV(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    N = cross(v2 - v1, v3 - v1, 2);
    nrm = sqrt(sum(N.^2,2));
    nrm(nrm==0) = 1;
    N = N ./ nrm;
end

function A = clippedAreaTriangleHalfSpace(p1,p2,p3, d1,d2,d3)
    pts = [p1; p2; p3];
    ds  = [d1; d2; d3];
    inside = ds <= 0;

    if all(~inside)
        A = 0; return
    elseif all(inside)
        A = triArea3D(p1,p2,p3); return
    end

    poly = [];
    for k = 1:3
        k2 = mod(k,3) + 1;
        Pk  = pts(k,:);  dk  = ds(k);
        Pk2 = pts(k2,:); dk2 = ds(k2);

        if dk <= 0
            poly = [poly; Pk];
        end

        if (dk <= 0 && dk2 > 0) || (dk > 0 && dk2 <= 0)
            t = dk / (dk - dk2);
            Pi = Pk + t*(Pk2 - Pk);
            poly = [poly; Pi];
        end
    end

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

function body = ensureBodyFields_fast(body)
    if ~isfield(body,'centroids')
        if isfield(body,'triCentroid')
            body.centroids = body.triCentroid;
        else
            error('Body missing centroids/triCentroid.');
        end
    end
    if ~isfield(body,'areas')
        if isfield(body,'triArea')
            body.areas = body.triArea;
        else
            error('Body missing areas/triArea.');
        end
    end
    if ~isfield(body,'bbox') || isempty(body.bbox)
        V = body.V;
        body.bbox = [min(V(:,1)) max(V(:,1)); min(V(:,2)) max(V(:,2)); min(V(:,3)) max(V(:,3))];
    end
end
