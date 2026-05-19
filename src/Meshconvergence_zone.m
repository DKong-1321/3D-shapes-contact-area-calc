% meshConvergence_3Pairs_localRefineSurface.m
% ------------------------------------------------------------
% Mesh convergence for 3 shape pairs using:
%   - Global coarse->fine surface generation via generateMesh (Hmax loop)
%   - Local surface refinement ("zoning in") ONLY on candidate triangles
%     near contact, via midpoint subdivision of selected triangles
%
% Outputs:
%   - Convergence curves vs mean edge length (NOT log faces)
%   - Saves Results to MAT + CSV in /out_meshConvergence_localRefine
%   - FINAL FIGURE per pair: mesh used for the final area calculation
%     with contact-region triangles highlighted.
%
% Notes:
%   - Local refinement here is surface-only and may be nonconforming at the
%     boundary (hanging nodes). For centroid/triangle-based contact-area
%     calculations this is usually fine.
% ------------------------------------------------------------

% ---- SAFE RENDERING (prevents WebGL/OpenGL crashes on heavy meshes) ----
set(groot, 'defaultFigureRenderer', 'painters');  % most stable
opengl('save','software');                        % forces software OpenGL

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Mesh convergence (3 pairs) with LOCAL surface refinement ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

%% ===== FILES (EDIT THESE EXACTLY TO YOUR FILENAMES) =====
cubeStlFine   = fullfile(modelDir,'50mm cube 4 iterations.stl');
sphereStlFine = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');
grooveStlFine = fullfile(modelDir,'Groove 22mm radius 4 it.stl');

assert(exist(cubeStlFine,'file')==2,   'Cube STL not found: %s', cubeStlFine);
assert(exist(sphereStlFine,'file')==2, 'Sphere STL not found: %s', sphereStlFine);
assert(exist(grooveStlFine,'file')==2, 'Groove STL not found: %s', grooveStlFine);

fprintf('STLs:\n');
fprintf('  Cube:   %s\n', cubeStlFine);
fprintf('  Sphere: %s\n', sphereStlFine);
fprintf('  Groove: %s\n\n', grooveStlFine);

%% ===== CONTACT SETTINGS =====
tol   = 0.20;         % mm contact tolerance used by computeContactArea_STS
pushIn = 0.5 * tol;   % small penetration for "overlap" placement

opts = struct();
opts.tol               = tol;
opts.roiExpandFactor   = 1.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 20;
opts.sampleMode        = 'centroid';
opts.sampleThreshold   = 0.5;

% If your computeContactArea_STS supports returning a mask, set this true.
% We still fall back to our candidate-mask highlighting if not available.
opts.returnMask        = true;

opts.useCentroidPrefilter = true;
opts.prefilterFactor      = 6.0;

% --- Debug / plotting (keep OFF for HPC) ---
opts.debug.enabled      = true;
opts.debug.everyLevel   = false;
opts.debug.levelsToPlot = [1 2 3];   % early levels = most informative
opts.debug.pause        = true;
opts.debug.showRoiBox   = true;

fprintf('Settings:\n');
fprintf('  tol=%.3f mm | pushIn=%.3f mm\n', tol, pushIn);
fprintf('  debug.enabled=%d\n\n', opts.debug.enabled);

%% ===== PLACEMENT SETTINGS =====
penCubeCube      = pushIn;
penSphereSphere  = 5.0;
penSphereGroove  = 24.0;
grooveSupportRxy = 60;      % mm support radius for groove top surface

%% ===== GLOBAL MESH CONVERGENCE (generateMesh Hmax) =====
HmaxStart  = 6.0;    % mm (coarse)
HmaxFactor = 0.70;   % multiply each step (smaller => finer)
HmaxMin    = 0.30;   % safety stop
maxLevels  = 25;     % safety cap

convTolPct = 5;      % convergence criterion on area change
convConsec = 1;      % consecutive hits needed
minLevels  = 5;      % don't stop too early

%% ===== LOCAL REFINEMENT CONTROLS =====
localRefIters = 3;      % local refinement cycles per global level
distFactor    = 1.0;    % candidate if dist <= distFactor * tol
normalDeg     = 60;     % normals: opposing-ish within this angle (deg)
% normalDeg larger => looser normal filter

%% ===== OUTPUTS =====
outDir = fullfile(projectRoot, 'out_meshConvergence_localRefine');
if ~exist(outDir,'dir'), mkdir(outDir); end
stamp  = datestr(now,'yyyymmdd_HHMMSS');
runTag = sprintf('meshConvLocal_tol%.3f_%s', tol, stamp);
matPath = fullfile(outDir, [runTag '.mat']);
csvPath = fullfile(outDir, [runTag '.csv']);

makeFigures = true;  % set false for no figure saving on HPC

%% ===== LOAD MASTERS (FIXED BODIES) =====
if ~exist('buildBodyStruct','file'), error('buildBodyStruct not found on path.'); end
if ~exist('computeContactArea_STS','file'), error('computeContactArea_STS not found on path.'); end

[FcF,VcF] = loadMeshAny(cubeStlFine);
[FsF,VsF] = loadMeshAny(sphereStlFine);
[FgF,VgF] = loadMeshAny(grooveStlFine);

% Keep your groove orientation
VgF = rotateAboutCentroid(VgF, rotxd(90));

cubeMaster   = buildBodyStruct(FcF,VcF);   cubeMaster.kdtree   = KDTreeSearcher(cubeMaster.centroids);
sphereMaster = buildBodyStruct(FsF,VsF);   sphereMaster.kdtree = KDTreeSearcher(sphereMaster.centroids);
grooveMaster = buildBodyStruct(FgF,VgF);   grooveMaster.kdtree = KDTreeSearcher(grooveMaster.centroids);

%% ===== RUN 3 PAIRS (LOCAL REFINEMENT) =====
Results = struct('name',{},'A',{},'Nfaces_moving',{},'dA_pct',{},'Hmax',{}, ...
                 'meanEdgeLen_mm',{},'medianEdgeLen_mm',{});

fprintf('--- Pair 1/3: cube-cube (local refine on moving cube) ---\n');
Results(1) = runPair_localRefineSurface( ...
    'cube-cube', cubeMaster, cubeStlFine, ...
    @(Vfix,Vmov) placeTwoBodies_overlap(Vfix,Vmov,'cube',penCubeCube), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, ...
    localRefIters, distFactor, normalDeg, ...
    convTolPct, convConsec, minLevels, tol, opts, 2000);

fprintf('\n--- Pair 2/3: sphere-sphere (local refine on moving sphere) ---\n');
Results(2) = runPair_localRefineSurface( ...
    'sphere-sphere', sphereMaster, sphereStlFine, ...
    @(Vfix,Vmov) placeTwoBodies_overlap(Vfix,Vmov,'sphere',penSphereSphere), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, ...
    localRefIters, distFactor, normalDeg, ...
    convTolPct, convConsec, minLevels, tol, opts, 2100);

fprintf('\n--- Pair 3/3: sphere-in-groove (local refine on moving sphere) ---\n');
Results(3) = runPair_localRefineSurface( ...
    'sphere-in-groove', grooveMaster, sphereStlFine, ...
    @(Vfix,Vmov) placeSphereInGroove_localSupport(Vfix,Vmov,penSphereGroove,grooveSupportRxy), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, ...
    localRefIters, distFactor, normalDeg, ...
    convTolPct, convConsec, minLevels, tol, opts, 2200);

%% ===== POST-PROCESS: TABLE + PLOTS + SAVE =====
T = buildResultsTable(Results);

writetable(T, csvPath);
save(matPath, 'Results', 'T', 'opts', 'tol', 'HmaxStart', 'HmaxFactor', 'HmaxMin', 'maxLevels', ...
    'localRefIters', 'distFactor', 'normalDeg', 'convTolPct', 'convConsec', 'minLevels', ...
    'cubeStlFine', 'sphereStlFine', 'grooveStlFine');

fprintf('\nSaved:\n  MAT: %s\n  CSV: %s\n', matPath, csvPath);

if makeFigures
    % Area vs mean edge
    f1 = figure('Color','w'); hold on; grid on
    xlabel('Mean surface edge length of moving mesh (mm)')
    ylabel('Contact area (mm^2)')
    title(sprintf('Mesh convergence (LOCAL refine), tol=%.3f mm', tol))
    for i = 1:numel(Results)
        plot(Results(i).meanEdgeLen_mm, Results(i).A, '-o', 'DisplayName', Results(i).name);
    end
    set(gca,'XDir','reverse'); legend('Location','best')
    saveas(f1, fullfile(outDir, [runTag '_A_vs_meanEdge.png']));

    % dA vs mean edge
    f2 = figure('Color','w'); hold on; grid on
    xlabel('Mean surface edge length of moving mesh (mm)')
    ylabel('% change vs previous level')
    title('Convergence rate (% change)')
    for i = 1:numel(Results)
        plot(Results(i).meanEdgeLen_mm, Results(i).dA_pct, '-o', 'DisplayName', Results(i).name);
    end
    set(gca,'XDir','reverse'); legend('Location','best')
    saveas(f2, fullfile(outDir, [runTag '_dA_vs_meanEdge.png']));
end

fprintf('\n=== Final convergence summary ===\n');
for i = 1:numel(Results)
    r = Results(i);
    fprintf('%-16s\n', r.name);
    fprintf('  Hmax:     %s\n', mat2str(round(r.Hmax,3)));
    fprintf('  faces:    %s\n', mat2str(r.Nfaces_moving));
    fprintf('  meanEdge: %s\n', mat2str(round(r.meanEdgeLen_mm,3)));
    fprintf('  area:     %s\n', mat2str(round(r.A,3)));
    fprintf('  dA%%:      %s\n\n', mat2str(round(r.dA_pct,3)));
end

%% ============================ FUNCTIONS ============================

function [F,V] = loadMeshAny(path)
    if exist('loadStlMesh','file')
        [F,V] = loadStlMesh(path);
    elseif exist('loadAnyStl','file')
        [F,V] = loadAnyStl(path);
    else
        error('Neither loadStlMesh nor loadAnyStl found on path.');
    end
end

function out = runPair_localRefineSurface(name, masterBody0, movStlPath, placer, ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, ...
    localRefIters, distFactor, normalDeg, ...
    convTolPct, convConsec, minLevels, tol, opts, figBase)

    % Preallocate per global level
    A   = nan(1,maxLevels);
    Nf  = nan(1,maxLevels);
    dA  = nan(1,maxLevels);
    Hs  = nan(1,maxLevels);
    Em  = nan(1,maxLevels);
    E50 = nan(1,maxLevels);

    Vfix0 = masterBody0.V;

    % Ensure master KD tree exists
    if ~isfield(masterBody0,'kdtree') || isempty(masterBody0.kdtree)
        masterBody0.kdtree = KDTreeSearcher(masterBody0.centroids);
    end

    cosThresh = cosd(normalDeg);

    streak = 0;
    H = HmaxStart;

    % Track final for the requested figure
    finalMaster = [];
    finalMov = [];
    finalMask = [];
    finalArea = NaN;

    for k = 1:maxLevels
        if H < HmaxMin
            fprintf('  Reached HmaxMin=%.3f mm, stopping.\n', HmaxMin);
            A   = A(1:max(0,k-1));
            Nf  = Nf(1:max(0,k-1));
            dA  = dA(1:max(0,k-1));
            Hs  = Hs(1:max(0,k-1));
            Em  = Em(1:max(0,k-1));
            E50 = E50(1:max(0,k-1));
            break
        end

        % ---- (1) Coarse surface from STL (global level) ----
        [Fmov, Vmov] = surfaceFromStl_generateMesh(movStlPath, H);

        % ---- (2) Placement ----
        [VfixP, VmovP] = placer(Vfix0, Vmov);

        masterBody = masterBody0;
        masterBody.V    = VfixP;
        masterBody.bbox = bboxFromV(VfixP);

        movBody = buildBodyStruct(Fmov, VmovP);
        movBody.bbox = bboxFromV(VmovP);

        % ---- (3) Local refinement loop on moving surface ----
        for it = 1:localRefIters
            candMask = computeCandidateTriangles(masterBody, movBody, tol, distFactor, cosThresh);

            if ~any(candMask)
                break
            end

            [Fref, Vref] = refineTrianglesMidpoint(movBody.F, movBody.V, candMask);
            movBody = buildBodyStruct(Fref, Vref);
            movBody.bbox = bboxFromV(Vref);

            if isfield(opts,'debug') && opts.debug.enabled
                doPlot = opts.debug.everyLevel || any(k == opts.debug.levelsToPlot);
                if doPlot
                    debugPlotPair(figBase + 1000 + 10*k + it, ...
                        sprintf('%s | global %d | local %d',name,k,it), it, ...
                        masterBody, movBody, tol, opts);
                    if isfield(opts.debug,'pause') && opts.debug.pause
                        waitforbuttonpress;
                    end
                end
            end
        end

        % ---- (4) Contact computation on locally-refined surface ----
        res = computeContactArea_STS(masterBody, movBody, tol, opts);

        contactMask = [];
        if isstruct(res)
            if isfield(res,'contactArea'), A(k) = res.contactArea; else, A(k) = NaN; end
            if isfield(res,'contactMask'), contactMask = res.contactMask; end
            if isempty(contactMask) && isfield(res,'maskSlave'), contactMask = res.maskSlave; end
            if isempty(contactMask) && isfield(res,'mask'), contactMask = res.mask; end
        else
            A(k) = res;
        end

        % Fallback: if algorithm doesn't return a mask, highlight candidates
        if isempty(contactMask)
            contactMask = computeCandidateTriangles(masterBody, movBody, tol, 1.0, cosThresh);
        end

        % ---- (5) Logging metrics ----
        Nf(k) = size(movBody.F,1);
        Hs(k) = H;
        [Em(k), E50(k)] = meshEdgeStats(movBody.F, movBody.V);

        if k >= 2 && ~isnan(A(k-1)) && A(k-1) ~= 0
            dA(k) = 100 * abs(A(k) - A(k-1)) / abs(A(k-1));
        end

        if k == 1
            fprintf('  Level %d | Hmax=%.3f | faces=%d | meanEdge=%.3f | A=%.4f\n', ...
                k, H, Nf(k), Em(k), A(k));
        else
            fprintf('  Level %d | Hmax=%.3f | faces=%d | meanEdge=%.3f | A=%.4f | dA%%=%.3f\n', ...
                k, H, Nf(k), Em(k), A(k), dA(k));
        end

        % Store final state each level (final one will be plotted)
        finalMaster = masterBody;
        finalMov    = movBody;
        finalMask   = contactMask;
        finalArea   = A(k);

        % ---- (6) Convergence ----
        if k >= 2 && ~isnan(dA(k))
            if k >= minLevels && dA(k) < convTolPct
                streak = streak + 1;
            else
                streak = 0;
            end
            if k >= minLevels && streak >= convConsec
                fprintf('  -> Converged: last %d dA%% < %.2f%%. Stopping.\n', convConsec, convTolPct);
                A   = A(1:k);
                Nf  = Nf(1:k);
                dA  = dA(1:k);
                Hs  = Hs(1:k);
                Em  = Em(1:k);
                E50 = E50(1:k);
                break
            end
        end

        H = H * HmaxFactor;
    end

    % ---- FINAL FIGURE requested: mesh used for final area calculation ----
    showFinalAreaMesh(figBase + 9000, name, finalMaster, finalMov, finalMask, tol, finalArea);

    out = struct('name',name,'A',A,'Nfaces_moving',Nf,'dA_pct',dA,'Hmax',Hs, ...
                 'meanEdgeLen_mm',Em,'medianEdgeLen_mm',E50);
end

function mask = computeCandidateTriangles(masterBody, movBody, tol, distFactor, cosThresh)
% Candidate triangles on moving surface based on:
%  - ROI overlap (expanded)
%  - KD distance to master triangle centroids
%  - Normal alignment (opposing-ish)

    C = movBody.centroids; % Nx3

    % ROI bbox overlap
    bbM = masterBody.bbox;
    bbS = movBody.bbox;
    roiMin = max(bbM(:,1), bbS(:,1)) - 2*tol;
    roiMax = min(bbM(:,2), bbS(:,2)) + 2*tol;

    inRoi = C(:,1) >= roiMin(1) & C(:,1) <= roiMax(1) & ...
            C(:,2) >= roiMin(2) & C(:,2) <= roiMax(2) & ...
            C(:,3) >= roiMin(3) & C(:,3) <= roiMax(3);

    if ~any(inRoi)
        mask = false(size(C,1),1);
        return
    end

    if ~isfield(masterBody,'kdtree') || isempty(masterBody.kdtree)
        masterBody.kdtree = KDTreeSearcher(masterBody.centroids);
    end

    [im, d] = knnsearch(masterBody.kdtree, C(inRoi,:), 'K', 1);

    near = false(size(C,1),1);
    near(inRoi) = (d <= distFactor * tol);

    % Normal alignment (prefer opposing normals)
    if isfield(masterBody,'normals') && isfield(movBody,'normals') && ...
       ~isempty(masterBody.normals) && ~isempty(movBody.normals)

        idx = find(near);
        if isempty(idx)
            mask = near;
            return
        end

        nm = masterBody.normals(im,:);       % master normals nearest each inRoi centroid
        nmFull = nan(size(C,1),3);
        nmFull(inRoi,:) = nm;               % place back

        nmSel = nmFull(idx,:);
        nsSel = movBody.normals(idx,:);

        % dot ~ -1 => facing. Use (-dot) >= cosThresh
        align = (-sum(nmSel.*nsSel,2)) >= cosThresh;

        mask = false(size(C,1),1);
        mask(idx) = align;
    else
        mask = near;
    end
end

function [F2, V2] = refineTrianglesMidpoint(F, V, refineMask)
% Refine ONLY triangles where refineMask==true by splitting into 4 via edge midpoints.
% Nonconforming refinement at borders (acceptable for triangle-based contact workflows).

    idxRef  = find(refineMask);
    idxKeep = find(~refineMask);

    Fkeep = F(idxKeep,:);
    Fr = F(idxRef,:);
    nR = size(Fr,1);

    V2 = V;
    Fnew = zeros(4*nR, 3);

    for i = 1:nR
        a = Fr(i,1); b = Fr(i,2); c = Fr(i,3);

        vab = 0.5*(V2(a,:) + V2(b,:));
        vbc = 0.5*(V2(b,:) + V2(c,:));
        vca = 0.5*(V2(c,:) + V2(a,:));

        iab = size(V2,1)+1; V2(iab,:) = vab;
        ibc = size(V2,1)+1; V2(ibc,:) = vbc;
        ica = size(V2,1)+1; V2(ica,:) = vca;

        base = 4*(i-1);
        Fnew(base+1,:) = [a, iab, ica];
        Fnew(base+2,:) = [iab, b, ibc];
        Fnew(base+3,:) = [ica, ibc, c];
        Fnew(base+4,:) = [iab, ibc, ica];
    end

    F2 = [Fkeep; Fnew];

    % Remove unused vertices
    [V2, F2] = removeUnusedVertices(V2, F2);
end

function [V2, F2] = removeUnusedVertices(V, F)
    used = false(size(V,1),1);
    used(F(:)) = true;
    map = zeros(size(V,1),1);
    map(used) = 1:nnz(used);
    V2 = V(used,:);
    F2 = map(F);
end

function [Fsurf, Vsurf] = surfaceFromStl_generateMesh(stlPath, Hmax)
% generateMesh creates a tet mesh from STL geometry; extract surface triangles via freeBoundary

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
        T = E(1:4,:).';        % take corner nodes of higher-order tets
    elseif nper == 3
        Fsurf = E.';           % already surface tris
        Vsurf = V;
        return
    else
        error('Unexpected msh.Elements size: %dx%d', size(E,1), size(E,2));
    end

    TR = triangulation(T, V);
    [Fsurf, Vsurf] = freeBoundary(TR);
end

function [meanE, medE] = meshEdgeStats(F, V)
% Mean and median unique edge length of a tri mesh

    if isempty(F) || isempty(V)
        meanE = NaN; medE = NaN; return
    end

    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');

    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2, 2));

    meanE = mean(L);
    medE  = median(L);
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
end

function [VfixP, VmovP] = placeTwoBodies_overlap(Vfix, Vmov, shapeType, penetration)
% Rough placement: align centroids in XY then set Z overlap by penetration

    VfixP = Vfix;

    cFix = mean(Vfix,1);
    cMov = mean(Vmov,1);
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

    cFix = mean(VfixP,1);
    cMov = mean(VmovP,1);

    dz = targetCenterDist - (cMov(3) - cFix(3));
    VmovP(:,3) = VmovP(:,3) + dz;
end

function [VfixP, VmovP] = placeSphereInGroove_localSupport(Vfix, Vmov, penetration, Rxy)
% Align XY centroids then sit the sphere onto the groove using a local Z support

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

function showFinalAreaMesh(figId, pairName, masterBody, movBody, contactMask, tol, A)
% Final display: mesh used for final area calc + contact region highlighted

    figure(figId); clf; set(gcf,'Color','w','Name',[pairName ' | FINAL AREA MESH']);
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | FINAL | tol=%.3f mm | A=%.4f mm^2', pairName, tol, A), 'Interpreter','none');

    patch('Faces', masterBody.F, 'Vertices', masterBody.V, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.12);

    patch('Faces', movBody.F, 'Vertices', movBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0.6 0], 'EdgeAlpha', 0.25);

    if ~isempty(contactMask) && any(contactMask)
        idx = find(contactMask);
        patch('Faces', movBody.F(idx,:), 'Vertices', movBody.V, ...
            'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end

    view(3); camlight; lighting gouraud
end

function debugPlotPair(figId, pairName, levelIdx, fixedBody, movBody, tol, opts)
    figure(figId); clf
    set(gcf,'Name',sprintf('%s | level %d', pairName, levelIdx));
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | level %d | tol=%.3f', pairName, levelIdx, tol), 'Interpreter','none');

    patch('Faces', fixedBody.F, 'Vertices', fixedBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.35);

    patch('Faces', movBody.F, 'Vertices', movBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0.7 0], 'EdgeAlpha', 0.35);

    view(3);
end

function T = buildResultsTable(Results)
% One row per (pair, level)

    rows = {};
    for i = 1:numel(Results)
        r = Results(i);
        n = numel(r.A);
        for k = 1:n
            rows(end+1,:) = {r.name, k, r.Hmax(k), r.Nfaces_moving(k), ...
                r.meanEdgeLen_mm(k), r.medianEdgeLen_mm(k), r.A(k), r.dA_pct(k)}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', ...
        {'pair','level','Hmax_mm','Nfaces','meanEdge_mm','medianEdge_mm','contactArea_mm2','dA_pct'});
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
