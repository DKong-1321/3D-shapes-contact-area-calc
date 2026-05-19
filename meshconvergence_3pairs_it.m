% meshConvergence_localRefine_contactEdgeOnly_showEveryLevel.m
% -------------------------------------------------------------------------
% Mesh convergence (3 pairs) using LOCAL refinement into contact patch:
%   - generateMesh ONCE to get an initial coarse MOVING surface
%   - compute PENETRATION contact area + get contact triangle mask (CLIPPED AREA)
%   - refine ONLY the contact triangles (midpoint subdivision)
%   - stop when dA% < convTolPct for convConsec consecutive levels (after minLevels)
%
% KEY UPDATE (your request):
%   Sphere–sphere placement is now:
%     - aligned in X/Y
%     - stacked along +Z (moving sphere starts above fixed sphere)
%     - then pushed down by "penetration" (mm)
%   => penetration=0 means just-touching; penetration>0 means overlap.
%
% CHANGES vs your current version:
%   (1) meanEdgeLen_mm computed ONLY on edges belonging to CONTACT triangles at THAT level
%       edges from Fmov(mask,:) only (not the whole surface).
%   (2) Shows a figure for EVERY mesh level (per pair), no clicking.
%
% Saves Results to:
%   <projectRoot>/results_meshConvergence_contactEdgeOnly.mat
%
% REQUIREMENTS on path (src/):
%   - buildBodyStruct
%   - loadStlMesh OR loadAnyStl
% -------------------------------------------------------------------------

clear; clc; close all;

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Mesh convergence (LOCAL refine into contact) | contact-only edge length | show EVERY level ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

%% ---- FILES (EDIT THESE EXACTLY TO YOUR FILENAMES) ----
cubeStlFine   = fullfile(modelDir,'50mm cube 4 iterations.stl');      % MASTER
sphereStlFine = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');   % MASTER
grooveStlFine = fullfile(modelDir,'hemisphere groove 50mm dia.stl');     % MASTER

assert(exist(cubeStlFine,'file')==2,   'Cube STL not found: %s', cubeStlFine);
assert(exist(sphereStlFine,'file')==2, 'Sphere STL not found: %s', sphereStlFine);
assert(exist(grooveStlFine,'file')==2, 'Groove STL not found: %s', grooveStlFine);

fprintf('STLs:\n');
fprintf('  Cube:   %s\n', cubeStlFine);
fprintf('  Sphere: %s\n', sphereStlFine);
fprintf('  Groove: %s\n\n', grooveStlFine);

%% ---- CONTACT + PLACEMENT SETTINGS ----
tol = 0.20; % mm (label / initial tol used during convergence run)

opts = struct();
opts.tol = tol;

% Debug plotting (positions) - NO CLICKS
opts.debugPlot = true;
opts.debugPlotEveryLevel = true;   % <-- SHOW EVERY LEVEL
opts.debugPlotPause = false;
opts.showRoiBox = false;

% Penetrations for placement (mm)
% NOTE: For two identical 50mm-diameter spheres, valid penetration range is [0..50].
% penetration=0 => just touching (stacked in +Z), penetration>0 => overlap.
penCubeCube     = 1;
penSphereSphere = 5;   % <-- set to what you actually want, e.g. 5 mm
penSphereGroove = 28;

grooveSupportRxy = 60;

fprintf('Settings:\n');
fprintf('  tol(label) = %.3f mm\n', tol);
fprintf('  penetrations (mm): cube=%.3f sphere=%.3f groove=%.3f\n\n', penCubeCube, penSphereSphere, penSphereGroove);

%% ---- LOCAL REFINEMENT / CONVERGENCE SETTINGS ----
HmaxStart   = 6.0;   % generateMesh Hmax ONLY for initial coarse MOVING surface (mm)
maxLevels   = 25;    % safety cap

convTolPct  = 5;     % convergence when dA% < 5
convConsec  = 3;     % 3 in a row
minLevels   = 3;

doDryRun = false;    % placement-only check

%% ---- LOAD MASTERS (fixed bodies) ----
[FcF,VcF] = loadMeshAny(cubeStlFine);
[FsF,VsF] = loadMeshAny(sphereStlFine);
[FgF,VgF] = loadMeshAny(grooveStlFine);

% Groove rotation: lay it down
VgF = rotateAboutCentroid(VgF, rotxd(90));

cubeMaster0   = buildBodyStruct(FcF,VcF);
sphereMaster0 = buildBodyStruct(FsF,VsF);
grooveMaster0 = buildBodyStruct(FgF,VgF);

%% ---- DRY RUN (PLACEMENT ONLY) ----
if doDryRun
    fprintf('=== DRY RUN (placement only) ===\n');

    % Generate coarse moving surfaces (same as main script)
    [FmovC, VmovC] = surfaceFromStl_generateMesh(cubeStlFine,   HmaxStart);
    [FmovS, VmovS] = surfaceFromStl_generateMesh(sphereStlFine, HmaxStart);

    % ---------------- cube-cube ----------------
    [VfixP,VmovP] = placeTwoBodies_overlap(VfixFromBody(cubeMaster0), VmovC, 'cube', penCubeCube);
    fixedBody = buildBodyStruct(cubeMaster0.F, VfixP); fixedBody.bbox = bboxFromV(VfixP);
    movBody   = buildBodyStruct(FmovC, VmovP);         movBody.bbox   = bboxFromV(VmovP);
    debugPlotPair(1001,'cube-cube (dry run)',1,fixedBody,movBody,tol,opts);

    % ---------------- sphere-sphere (STACKED +Z) ----------------
    [VfixP,VmovP] = placeTwoSpheres_stackZ(VfixFromBody(sphereMaster0), VmovS, penSphereSphere);

    % ---- DIAGNOSTIC: measured penetration (mm) ----
    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(VmovP);
    Rfix = 0.5 * (bbFix(1,2) - bbFix(1,1));
    Rmov = 0.5 * (bbMov(1,2) - bbMov(1,1));
    cFix = bboxCenter(VfixP);
    cMov = bboxCenter(VmovP);
    penetration_measured_mm = (Rfix + Rmov) - norm(cMov - cFix);
    fprintf('\nSphere–sphere diagnostic:\n');
    fprintf('  Requested penetration = %.3f mm\n', penSphereSphere);
    fprintf('  Measured penetration  = %.6f mm\n\n', penetration_measured_mm);
    % -----------------------------------------------------------

    fixedBody = buildBodyStruct(sphereMaster0.F, VfixP); fixedBody.bbox = bboxFromV(VfixP);
    movBody   = buildBodyStruct(FmovS, VmovP);           movBody.bbox   = bboxFromV(VmovP);
    debugPlotPair(1002,'sphere-sphere (dry run)',1,fixedBody,movBody,tol,opts);

    % ---------------- sphere-in-groove ----------------
    [VfixP,VmovP] = placeSphereInGroove_localSupport(VfixFromBody(grooveMaster0), VmovS, penSphereGroove, grooveSupportRxy);
    fixedBody = buildBodyStruct(grooveMaster0.F, VfixP); fixedBody.bbox = bboxFromV(VfixP);
    movBody   = buildBodyStruct(FmovS, VmovP);           movBody.bbox   = bboxFromV(VmovP);
    debugPlotPair(1003,'sphere-in-groove (dry run)',1,fixedBody,movBody,tol,opts);

    fprintf('\nDry run done.\nSet doDryRun=false and rerun.\n');
    return
end

%% ---- RUN 3 PAIRS ----
Results = struct('name',{},'A',{},'Nfaces_moving',{},'dA_pct',{}, ...
                 'meanEdgeLen_contact_mm',{},'medianEdgeLen_contact_mm',{}, ...
                 'meanEdgeLen_all_mm',{}, ...
                 'F_final',{},'V_finalPlaced',{},'faceLvl_final',{}, ...
                 'F_master',{},'V_masterPlaced',{},'mask_last',{});

fprintf('--- Pair 1/3: cube-cube ---\n');
Results(1) = runPair_refineContactTris_only( ...
    'cube-cube', cubeMaster0, cubeStlFine, ...
    @(Vfix,Vmov) placeTwoBodies_overlap(Vfix,Vmov,'cube',penCubeCube), ...
    HmaxStart, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2000);

fprintf('\n--- Pair 2/3: sphere-sphere (STACKED +Z) ---\n');
Results(2) = runPair_refineContactTris_only( ...
    'sphere-sphere', sphereMaster0, sphereStlFine, ...
    @(Vfix,Vmov) placeTwoSpheres_stackZ(Vfix,Vmov,penSphereSphere), ...
    HmaxStart, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2100);

fprintf('\n--- Pair 3/3: sphere-in-groove ---\n');
Results(3) = runPair_refineContactTris_only( ...
    'sphere-in-groove', grooveMaster0, sphereStlFine, ...
    @(Vfix,Vmov) placeSphereInGroove_localSupport(Vfix,Vmov,penSphereGroove,grooveSupportRxy), ...
    HmaxStart, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2200);
%% ---- EXPORT FINAL SPHERE-GROOVE CONTACT FIGURES ----
% (1) groove transparent + sphere + contact highlighted
% (2) sphere only + contact highlighted

outFigDir = fullfile(projectRoot,'results_figs');
if ~exist(outFigDir,'dir'), mkdir(outFigDir); end

R = Results(3); % sphere-in-groove

% Ensure we have a usable mask
mask = [];
if isfield(R,'mask_last'), mask = R.mask_last; end
if isempty(mask) || ~islogical(mask) || numel(mask) ~= size(R.F_final,1)
    warning('sphere-in-groove: mask_last missing/invalid; figures will be exported without contact highlight.');
    mask = false(size(R.F_final,1),1);
end

% ---------- Figure 1: Groove transparent + Sphere + Contact ----------
fig1 = figure('Color','w','Name','sphere-in-groove | groove transparent + contact');
hold on; axis equal; grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('sphere-in-groove | FINAL | tol(label)=%.3f mm', tol), 'Interpreter','none');

% Groove (master): transparent surface, no mesh edges
patch('Faces', R.F_master, 'Vertices', R.V_masterPlaced, ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.08, ...
    'EdgeColor', 'none');

% Sphere (moving): show mesh lightly
patch('Faces', R.F_final, 'Vertices', R.V_finalPlaced, ...
    'FaceColor', [0 0.7 0], 'FaceAlpha', 0.10, ...
    'EdgeColor', [0 0 0], 'EdgeAlpha', 0.18, ...
    'LineWidth', 0.25);

% Contact triangles (highlight): solid magenta
if any(mask)
    patch('Faces', R.F_final(mask,:), 'Vertices', R.V_finalPlaced, ...
        'FaceColor', [1 0 1], 'FaceAlpha', 0.95, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, ...
        'LineWidth', 0.30);
end

view(3); camlight; lighting gouraud

fpath1_png = fullfile(outFigDir,'sphere_in_groove_FINAL_grooveTransparent_contact.png');
fpath1_fig = fullfile(outFigDir,'sphere_in_groove_FINAL_grooveTransparent_contact.fig');
exportgraphics(fig1, fpath1_png, 'Resolution', 300);
savefig(fig1, fpath1_fig);

% ---------- Figure 2: Sphere only + Contact ----------
fig2 = figure('Color','w','Name','sphere-in-groove | sphere only + contact');
hold on; axis equal; grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('sphere only | FINAL contact | tol(label)=%.3f mm', tol), 'Interpreter','none');

% Sphere only: show mesh lightly
patch('Faces', R.F_final, 'Vertices', R.V_finalPlaced, ...
    'FaceColor', [0 0.7 0], 'FaceAlpha', 0.10, ...
    'EdgeColor', [0 0 0], 'EdgeAlpha', 0.18, ...
    'LineWidth', 0.25);

% Contact highlight
if any(mask)
    patch('Faces', R.F_final(mask,:), 'Vertices', R.V_finalPlaced, ...
        'FaceColor', [1 0 1], 'FaceAlpha', 0.95, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, ...
        'LineWidth', 0.30);
end

view(3); camlight; lighting gouraud

fpath2_png = fullfile(outFigDir,'sphere_only_FINAL_contact.png');
fpath2_fig = fullfile(outFigDir,'sphere_only_FINAL_contact.fig');
exportgraphics(fig2, fpath2_png, 'Resolution', 300);
savefig(fig2, fpath2_fig);

fprintf('\nExported final contact figures to:\n  %s\n  %s\n  %s\n  %s\n\n', ...
    fpath1_png, fpath1_fig, fpath2_png, fpath2_fig);

%% ---- COMBINED PLOTS (LOG CONTACT-PATCH EDGE LENGTH) ----
figure('Color','w'); hold on; grid on
xlabel('log_{10}(mean edge length of CONTACT triangles)  [mm]');
ylabel('Penetration contact area (mm^2)')
title(sprintf('Mesh convergence (refine contact tris only) | tol(label)=%.3f mm', tol))
for i = 1:numel(Results)
    x = log10(Results(i).meanEdgeLen_contact_mm);
    y = Results(i).A;
    plot(x, y, '-o', 'DisplayName', Results(i).name);
end
set(gca,'XDir','reverse')
legend('Location','best')

figure('Color','w'); hold on; grid on
xlabel('log_{10}(mean edge length of CONTACT triangles)  [mm]');
ylabel('% change vs previous level')
title('Convergence rate (% change)')
for i = 1:numel(Results)
    x = log10(Results(i).meanEdgeLen_contact_mm);
    y = Results(i).dA_pct;
    plot(x, y, '-o', 'DisplayName', Results(i).name);
end
set(gca,'XDir','reverse')
legend('Location','best')

%% ---- FINAL ZONING FIGURES ----
for i = 1:numel(Results)
    plotFinalZoningPair(3000+i, Results(i), tol);
end

fprintf('\n=== Final convergence summary ===\n');
for i = 1:numel(Results)
    r = Results(i);
    fprintf('%-16s | levels=%d\n', r.name, numel(r.A));
    fprintf('  faces: %s\n', mat2str(r.Nfaces_moving));
    fprintf('  mean edge CONTACT (mm): %s\n', mat2str(round(r.meanEdgeLen_contact_mm,3)));
    fprintf('  areas: %s\n', mat2str(round(r.A,3)));
    fprintf('  dA%% :  %s\n\n', mat2str(round(r.dA_pct,3)));
end
fprintf('=== DONE ===\n');

%% ---- SAVE converged meshes + logs for tol sweep script ----
outPath = fullfile(projectRoot,'results_meshConvergence_contactEdgeOnly.mat');
save(outPath,'Results','tol','HmaxStart','convTolPct','convConsec','minLevels','maxLevels');
fprintf('Saved: %s\n', outPath);

%% ============================ FUNCTIONS ============================

function Vfix = VfixFromBody(body)
    Vfix = body.V;
end

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

    % CONTACT-PATCH edge length (what you want for log plot)
    EmC  = nan(1,maxLevels);
    E50C = nan(1,maxLevels);

    % Optional: whole mesh edge length (kept as a sanity reference)
    EmAll = nan(1,maxLevels);

    % --- initial coarse moving surface (UNPLACED coordinates) ---
    [Fmov, Vmov] = surfaceFromStl_generateMesh(movStlPath, HmaxStart);
    faceLvl = zeros(size(Fmov,1),1,'uint16');

    % PLACE ONCE, FREEZE TRANSFORM
    [VfixP0, VmovP0] = placer(Vfix0, Vmov);

    % Frozen translation for moving mesh (use bbox-centre to avoid drift)
    tMove = bboxCenter(VmovP0) - bboxCenter(Vmov);

    % Frozen master placement
    VfixPlacedConst = VfixP0;

    % Store placed versions
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

        % Debug plot EVERY level
        if isfield(opts,'debugPlot') && opts.debugPlot
            if isfield(opts,'debugPlotEveryLevel') && opts.debugPlotEveryLevel
                debugPlotPair(figBase + k, name, k, masterBody, movBody, tol, opts);
            end
        end

        % Contact / penetration area + mask
        res = computeContactArea_STS_penetration(masterBody, movBody, tol, opts);
        mask = res.contactMask;
        A(k) = res.contactArea;

        if isempty(mask) || ~islogical(mask) || numel(mask) ~= size(Fmov,1)
            mask = false(size(Fmov,1),1);
            warning('%s: contact mask not returned/invalid at level %d.', name, k);
        end
        mask_last = mask;

        % CONTACT-only edge length (key)
        if any(mask)
            [EmC(k), E50C(k)] = meshEdgeStats_faces(Fmov(mask,:), VmovP);
        else
            EmC(k)  = NaN;
            E50C(k) = NaN;
        end

        % Whole-mesh mean edge (optional sanity)
        EmAll(k) = meshEdgeStats_all(Fmov, VmovP);

        % Convergence metric on area
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

    % Show COARSE / STEP BEFORE FINAL / FINAL
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
        Fsurf = E.'; Vsurf = V;
        return
    else
        error('Unexpected msh.Elements size: %dx%d', size(E,1), size(E,2));
    end

    TR = triangulation(T, V);
    [Fsurf, Vsurf] = freeBoundary(TR);
end

function meanE = meshEdgeStats_all(F, V)
    if isempty(F) || isempty(V)
        meanE = NaN; return
    end
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2,2));
    meanE = mean(L);
end

function [meanE, medE] = meshEdgeStats_faces(Fsub, V)
    if isempty(Fsub) || isempty(V)
        meanE = NaN; medE = NaN; return
    end
    E = [Fsub(:,[1 2]); Fsub(:,[2 3]); Fsub(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2,2));
    meanE = mean(L);
    medE  = median(L);
end

function [VfixP, VmovP] = placeTwoBodies_overlap(Vfix, Vmov, shapeType, penetration)
% placeTwoBodies_overlap
% Generic block/sphere placer (kept for cube-cube).
% For spheres, this aligns centroids and sets centre distance in Z.
% (You are NOT using this for sphere-sphere anymore: use placeTwoSpheres_stackZ)

    VfixP = Vfix;

    cFix = mean(Vfix, 1);
    cMov = mean(Vmov, 1);
    VmovP = Vmov + (cFix - cMov);

    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(VmovP);

    if strcmpi(shapeType,'sphere')
        Rfix = 0.5 * (bbFix(1,2) - bbFix(1,1));
        Rmov = 0.5 * (bbMov(1,2) - bbMov(1,1));
        targetCenterDist = (Rfix + Rmov) - penetration;
    else
        Hfix = 0.5 * (bbFix(3,2) - bbFix(3,1));
        Hmov = 0.5 * (bbMov(3,2) - bbMov(3,1));
        targetCenterDist = (Hfix + Hmov) - penetration;
    end

    cFix = mean(VfixP, 1);
    cMov = mean(VmovP, 1);

    currentDz = (cMov(3) - cFix(3));
    dz = targetCenterDist - currentDz;

    VmovP(:,3) = VmovP(:,3) + dz;
end

function [VfixP, VmovP] = placeTwoSpheres_stackZ(Vfix, Vmov, penetration)
% placeTwoSpheres_stackZ
% - Align X/Y using bbox centres
% - Stack moving sphere ABOVE fixed sphere along +Z
% - Apply penetration by reducing centre distance by "penetration"
%
% penetration=0  -> just touching
% penetration>0  -> overlap
% penetration max (equal spheres) -> Rfix+Rmov

    VfixP = Vfix;

    cFix = bboxCenter(VfixP);
    cMov = bboxCenter(Vmov);

    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(Vmov);
    Rfix = 0.5 * (bbFix(1,2) - bbFix(1,1));
    Rmov = 0.5 * (bbMov(1,2) - bbMov(1,1));

    penetration = max(0, min(penetration, (Rfix + Rmov)));

    % align XY
    VmovP = Vmov;
    VmovP(:,1) = VmovP(:,1) + (cFix(1) - cMov(1));
    VmovP(:,2) = VmovP(:,2) + (cFix(2) - cMov(2));

    % recompute centre
    cMov2 = bboxCenter(VmovP);

    targetCenterDist = (Rfix + Rmov) - penetration;

    % stack above (+Z)
    dz = (cFix(3) + targetCenterDist) - cMov2(3);
    VmovP(:,3) = VmovP(:,3) + dz;
end

function [VfixP, VmovP] = placeSphereInGroove_localSupport(Vfix, Vmov, penetration, Rxy)
    VfixP = Vfix;

    cFixXY = mean(Vfix(:,1:2),1);
    cMovXY = mean(Vmov(:,1:2),1);
    VmovP = Vmov;
    VmovP(:,1:2) = VmovP(:,1:2) + (cFixXY - cMovXY);

    cXY = mean(VmovP(:,1:2),1);
    dxy = hypot(VfixP(:,1) - cXY(1), VfixP(:,2) - cXY(2));
    idx = dxy <= Rxy;

    if any(idx)
        zSupport = max(VfixP(idx,3));
    else
        zSupport = max(VfixP(:,3));
    end

    zMovBottom = min(VmovP(:,3));
    VmovP(:,3) = VmovP(:,3) + ((zSupport - penetration) - zMovBottom);
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
