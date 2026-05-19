% meshConvergence_3Pairs_generateMesh.m
% Mesh convergence (3 pairs) using generateMesh (tet mesh -> surface freeBoundary)
% Hygiene/perf updates:
%  - Debug plotting OFF by default (HPC-friendly)
%  - Debug plotting controlled via opts.debug
%  - Pipeline split: mesh -> placement -> contact -> logging -> postprocess
%  - Results saved to MAT + CSV automatically, figures optional

clear; clc;

%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Mesh convergence (3 pairs): generateMesh coarse->fine until converged ===\n');
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

%% ===== CONTACT + PLACEMENT SETTINGS =====
tol = 0.20;                 % mm
pushIn = 0.5 * tol;

opts = struct();
opts.tol               = tol;
opts.roiExpandFactor   = 1.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 20;
opts.sampleMode        = 'centroid';
opts.sampleThreshold   = 0.5;
opts.returnMask        = false;
opts.useCentroidPrefilter = true;
opts.prefilterFactor      = 6.0;

% --- DEBUG SETTINGS (HPC-friendly defaults) ---
opts.debug = struct();
opts.debug.enabled     = false;   % <<< set true only when you want figures
opts.debug.everyLevel  = false;   % if true, plot every level (slow)
opts.debug.levelsToPlot = 1;      % e.g. [1 3 5] plots those levels only
opts.debug.pause       = false;   % if true, waits for keypress
opts.debug.showRoiBox  = false;

penCubeCube     = pushIn;
penSphereSphere = 5.0;
penSphereGroove = 24.0;
grooveSupportRxy = 60;      % mm

fprintf('Settings:\n');
fprintf('  tol = %.3f mm\n', tol);
fprintf('  pushIn = %.3f mm\n', pushIn);
fprintf('  debug.enabled=%d (everyLevel=%d pause=%d)\n\n', ...
    opts.debug.enabled, opts.debug.everyLevel, opts.debug.pause);

%% ===== MESH CONVERGENCE SETTINGS (GENERATEMESH) =====
HmaxStart   = 6.0;      % mm
HmaxFactor  = 0.70;     % refinement multiplier
HmaxMin     = 0.30;     % safety stop
maxLevels   = 25;       % safety cap

convTolPct  = 5;        % convergence criterion
convConsec  = 1;        % consecutive hits required (set 3 if you really want 3)
minLevels   = 5;        % don't stop too early

doDryRun = false;       % placement-only mode

%% ===== OUTPUT SETTINGS (AUTO-SAVE) =====
outDir = fullfile(projectRoot, 'out_meshConvergence');
if ~exist(outDir,'dir'), mkdir(outDir); end

stamp = datestr(now,'yyyymmdd_HHMMSS');
runTag = sprintf('meshConv_tol%.3f_%s', tol, stamp);
matPath = fullfile(outDir, [runTag '.mat']);
csvPath = fullfile(outDir, [runTag '.csv']);

makeFigures = true;  % set false on HPC if you want zero figures saved

%% ===== LOAD MASTERS =====
[FcF,VcF] = loadMeshAny(cubeStlFine);
[FsF,VsF] = loadMeshAny(sphereStlFine);
[FgF,VgF] = loadMeshAny(grooveStlFine);

VgF = rotateAboutCentroid(VgF, rotxd(90));  % keep groove rotation

if ~exist('buildBodyStruct','file'), error('buildBodyStruct not found on path.'); end
if ~exist('computeContactArea_STS','file'), error('computeContactArea_STS not found on path.'); end

cubeMaster   = buildBodyStruct(FcF,VcF); cubeMaster.kdtree = KDTreeSearcher(cubeMaster.centroids);
sphereMaster = buildBodyStruct(FsF,VsF); sphereMaster.kdtree = KDTreeSearcher(sphereMaster.centroids);
grooveMaster = buildBodyStruct(FgF,VgF); grooveMaster.kdtree = KDTreeSearcher(grooveMaster.centroids);

%% ===== DRY RUN (PLACEMENT ONLY) =====
if doDryRun
    fprintf('=== DRY RUN (placement only) ===\n');
    [FmovC, VmovC] = surfaceFromStl_generateMesh(cubeStlFine, HmaxStart);
    [FmovS, VmovS] = surfaceFromStl_generateMesh(sphereStlFine, HmaxStart);

    [VfixP,VmovP] = placeTwoBodies_overlap(cubeMaster.V, VmovC, 'cube', penCubeCube);
    debugPlotPair(1001,'cube-cube (dry run)',1,buildBodyStruct(FcF,VfixP),buildBodyStruct(FmovC,VmovP),tol,opts);

    [VfixP,VmovP] = placeTwoBodies_overlap(sphereMaster.V, VmovS, 'sphere', penSphereSphere);
    debugPlotPair(1002,'sphere-sphere (dry run)',1,buildBodyStruct(FsF,VfixP),buildBodyStruct(FmovS,VmovP),tol,opts);

    [VfixP,VmovP] = placeSphereInGroove_localSupport(grooveMaster.V, VmovS, penSphereGroove, grooveSupportRxy);
    debugPlotPair(1003,'sphere-in-groove (dry run)',1,buildBodyStruct(FgF,VfixP),buildBodyStruct(FmovS,VmovP),tol,opts);

    fprintf('\nDry run done.\nSet doDryRun=false and rerun.\n');
    return
end

%% ===== RUN 3 PAIRS =====
Results = struct('name',{},'A',{},'Nfaces_moving',{},'dA_pct',{},'Hmax',{}, ...
                 'meanEdgeLen_mm',{},'medianEdgeLen_mm',{});

fprintf('--- Pair 1/3: cube-cube ---\n');
Results(1) = runPair_generateMeshUntilConverged_clean( ...
    'cube-cube', cubeMaster, cubeStlFine, ...
    @(Vfix,Vmov) placeTwoBodies_overlap(Vfix,Vmov,'cube',penCubeCube), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2000);

fprintf('\n--- Pair 2/3: sphere-sphere ---\n');
Results(2) = runPair_generateMeshUntilConverged_clean( ...
    'sphere-sphere', sphereMaster, sphereStlFine, ...
    @(Vfix,Vmov) placeTwoBodies_overlap(Vfix,Vmov,'sphere',penSphereSphere), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2100);

fprintf('\n--- Pair 3/3: sphere-in-groove ---\n');
Results(3) = runPair_generateMeshUntilConverged_clean( ...
    'sphere-in-groove', grooveMaster, sphereStlFine, ...
    @(Vfix,Vmov) placeSphereInGroove_localSupport(Vfix,Vmov,penSphereGroove,grooveSupportRxy), ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, convTolPct, convConsec, minLevels, tol, opts, 2200);

%% ===== POST-PROCESS + SAVE =====
[T, figHandles] = postprocess_and_save(Results, tol, outDir, runTag, makeFigures);

save(matPath, 'Results', 'opts', 'tol', 'HmaxStart', 'HmaxFactor', 'HmaxMin', 'maxLevels', ...
    'convTolPct', 'convConsec', 'minLevels', 'cubeStlFine', 'sphereStlFine', 'grooveStlFine', 'T');

writetable(T, csvPath);

fprintf('\nSaved:\n  MAT: %s\n  CSV: %s\n', matPath, csvPath);

fprintf('\n=== Final convergence summary ===\n');
for i = 1:numel(Results)
    r = Results(i);
    fprintf('%-16s\n', r.name);
    fprintf('  Hmax:       %s\n', mat2str(round(r.Hmax,3)));
    fprintf('  faces:      %s\n', mat2str(r.Nfaces_moving));
    fprintf('  meanEdge:   %s\n', mat2str(round(r.meanEdgeLen_mm,3)));
    fprintf('  areas:      %s\n', mat2str(round(r.A,3)));
    fprintf('  dA%%:        %s\n\n', mat2str(round(r.dA_pct,3)));
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

function out = runPair_generateMeshUntilConverged_clean(name, masterBody0, movStlPath, placer, ...
    HmaxStart, HmaxFactor, HmaxMin, maxLevels, convTolPct, convConsec, minLevels, tol, opts, figBase)

    % --- PREALLOCATE ---
    A   = nan(1,maxLevels);
    Nf  = nan(1,maxLevels);
    dA  = nan(1,maxLevels);
    Hs  = nan(1,maxLevels);
    Em  = nan(1,maxLevels);
    E50 = nan(1,maxLevels);

    Vfix0 = masterBody0.V;

    streak = 0;
    H = HmaxStart;

    for k = 1:maxLevels
        if H < HmaxMin
            fprintf('  Reached HmaxMin=%.3f mm, stopping.\n', HmaxMin);
            [A,Nf,dA,Hs,Em,E50] = trimToKminus1(A,Nf,dA,Hs,Em,E50,k);
            break
        end

        % ---- (1) MESH GENERATION ----
        [Fmov, Vmov] = surfaceFromStl_generateMesh(movStlPath, H);
        Nf(k) = size(Fmov,1);
        Hs(k) = H;
        [Em(k), E50(k)] = meshEdgeStats(Fmov, Vmov);

        fprintf('  Level %d | Hmax=%.3f | faces=%d | meanEdge=%.3f mm ... ', ...
            k, H, Nf(k), Em(k));

        % ---- (2) PLACEMENT ----
        [VfixP, VmovP] = placer(Vfix0, Vmov);

        masterBody = masterBody0;
        masterBody.V    = VfixP;
        masterBody.bbox = bboxFromV(VfixP);

        movBody = buildBodyStruct(Fmov, VmovP);
        movBody.bbox = bboxFromV(VmovP);

        % ---- (3) OPTIONAL DEBUG PLOT ----
        if isfield(opts,'debug') && opts.debug.enabled
            doPlot = opts.debug.everyLevel || any(k == opts.debug.levelsToPlot);
            if doPlot
                debugPlotPair(figBase + k, name, k, masterBody, movBody, tol, opts);
                if isfield(opts.debug,'pause') && opts.debug.pause
                    fprintf('Press any key in figure to continue...\n');
                    waitforbuttonpress;
                end
            end
        end

        % ---- (4) CONTACT COMPUTATION ----
        res = computeContactArea_STS(masterBody, movBody, tol, opts);
        if isstruct(res), A(k) = res.contactArea; else, A(k) = res; end

        % ---- (5) CONVERGENCE LOGGING ----
        if k >= 2 && A(k-1) ~= 0
            dA(k) = 100 * abs(A(k) - A(k-1)) / abs(A(k-1));
        end

        if k == 1
            fprintf('A=%.4f\n', A(k));
        else
            fprintf('A=%.4f | dA%%=%.3f\n', A(k), dA(k));
        end

        if k >= 2 && ~isnan(dA(k))
            if k >= minLevels && dA(k) < convTolPct
                streak = streak + 1;
            else
                streak = 0;
            end

            if k >= minLevels && streak >= convConsec
                fprintf('  -> Converged: last %d dA%% < %.2f%%. Stopping.\n', convConsec, convTolPct);
                [A,Nf,dA,Hs,Em,E50] = trimToK(A,Nf,dA,Hs,Em,E50,k);
                break
            end
        end

        H = H * HmaxFactor;
    end

    out = struct('name',name,'A',A,'Nfaces_moving',Nf,'dA_pct',dA,'Hmax',Hs, ...
                 'meanEdgeLen_mm',Em,'medianEdgeLen_mm',E50);
end

function [A,Nf,dA,Hs,Em,E50] = trimToK(A,Nf,dA,Hs,Em,E50,k)
    A=A(1:k); Nf=Nf(1:k); dA=dA(1:k); Hs=Hs(1:k); Em=Em(1:k); E50=E50(1:k);
end
function [A,Nf,dA,Hs,Em,E50] = trimToKminus1(A,Nf,dA,Hs,Em,E50,k)
    if k<=1
        A=[]; Nf=[]; dA=[]; Hs=[]; Em=[]; E50=[];
    else
        A=A(1:k-1); Nf=Nf(1:k-1); dA=dA(1:k-1); Hs=Hs(1:k-1); Em=Em(1:k-1); E50=E50(1:k-1);
    end
end

function [T, figHandles] = postprocess_and_save(Results, tol, outDir, runTag, makeFigures)
    % Build a tidy table for CSV + easier later analysis
    rows = {};
    for i = 1:numel(Results)
        r = Results(i);
        n = numel(r.A);
        for k = 1:n
            rows(end+1,:) = {r.name, k, r.Hmax(k), r.Nfaces_moving(k), r.meanEdgeLen_mm(k), ...
                             r.medianEdgeLen_mm(k), r.A(k), r.dA_pct(k)}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', ...
        {'pair','level','Hmax_mm','Nfaces','meanEdge_mm','medianEdge_mm','contactArea_mm2','dA_pct'});

    figHandles = [];

    if makeFigures
        % Area vs mean edge
        f1 = figure('Color','w'); hold on; grid on
        xlabel('Mean surface edge length of moving mesh (mm)')
        ylabel('Contact area (mm^2)')
        title(sprintf('Mesh convergence (generateMesh), tol=%.3f mm', tol))
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

        figHandles = [f1 f2];
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

    V = msh.Nodes.';   % Nx3
    E = msh.Elements;  % nNodesPerElem x nElems
    nper = size(E,1);

    if nper == 4
        T = E.';           % linear tets
    elseif nper > 4
        T = E(1:4,:).';    % take corners for quadratic tets
    elseif nper == 3
        Fsurf = E.'; Vsurf = V;
        return
    else
        error('Unexpected msh.Elements size: %dx%d', size(E,1), size(E,2));
    end

    TR = triangulation(T, V);
    [Fsurf, Vsurf] = freeBoundary(TR);
end

function [meanE, medE] = meshEdgeStats(F, V)
    if isempty(F) || isempty(V)
        meanE = NaN; medE = NaN; return
    end
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E, 2);
    E = unique(E, 'rows');
    d = V(E(:,1),:) - V(E(:,2),:);
    L = sqrt(sum(d.^2, 2));
    meanE = mean(L);
    medE  = median(L);
end

function [VfixP, VmovP] = placeTwoBodies_overlap(Vfix, Vmov, shapeType, penetration)
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

    if isfield(opts,'debug') && isfield(opts.debug,'showRoiBox') && opts.debug.showRoiBox
        roiMin = max(movBody.bbox(:,1), fixedBody.bbox(:,1)) - opts.roiExpandFactor * tol;
        roiMax = min(movBody.bbox(:,2), fixedBody.bbox(:,2)) + opts.roiExpandFactor * tol;
        drawBBox3D(roiMin, roiMax);
    end

    view(3);
end

function drawBBox3D(roiMin, roiMax)
    xmin = roiMin(1); ymin = roiMin(2); zmin = roiMin(3);
    xmax = roiMax(1); ymax = roiMax(2); zmax = roiMax(3);

    P = [ xmin ymin zmin;
          xmax ymin zmin;
          xmax ymax zmin;
          xmin ymax zmin;
          xmin ymin zmax;
          xmax ymin zmax;
          xmax ymax zmax;
          xmin ymax zmax ];

    E = [1 2; 2 3; 3 4; 4 1; ...
         5 6; 6 7; 7 8; 8 5; ...
         1 5; 2 6; 3 7; 4 8];

    for e = 1:size(E,1)
        p1 = P(E(e,1),:);
        p2 = P(E(e,2),:);
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'r-', 'LineWidth', 1.5);
    end
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
