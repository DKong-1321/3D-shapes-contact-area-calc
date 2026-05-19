% meshConvergence_3Pairs_generateMesh_full.m
clear; clc; close all

scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
addpath(genpath(fullfile(projectRoot,'src')))
modelDir = fullfile(projectRoot,'model');

fprintf('\n=== Mesh convergence (3 pairs) using generateMesh remeshing ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

cubeStl   = fullfile(modelDir,'cube_50mm.stl');
sphereStl = fullfile(modelDir,'sphere_50mm.stl');
grooveStl = fullfile(modelDir,'Groove 22mm radius 4 it.stl');

assert(exist(cubeStl,'file')==2,   'Cube STL not found: %s', cubeStl);
assert(exist(sphereStl,'file')==2, 'Sphere STL not found: %s', sphereStl);
assert(exist(grooveStl,'file')==2, 'Groove STL not found: %s', grooveStl);

fprintf('Found STL files:\n');
fprintf('  Cube:   %s\n', cubeStl);
fprintf('  Sphere: %s\n', sphereStl);
fprintf('  Groove: %s\n\n', grooveStl);

tol = 0.20;   % mm

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

opts.debugPlot = true;
opts.debugPlotEveryLevel = false;
opts.debugPlotPause = true;
opts.plotCentroids = false;

fprintf('Settings:\n');
fprintf('  tol = %.3f mm\n', tol);
fprintf('  maxNeighbours = %d\n', opts.maxNeighbours);
fprintf('  ROI expand factor = %.2f\n', opts.roiExpandFactor);
fprintf('  centroid prefilter = %d (gate=%.2f*tol)\n', opts.useCentroidPrefilter, opts.prefilterFactor);
fprintf('  debugPlot = %d (everyLevel=%d, pause=%d)\n\n', opts.debugPlot, opts.debugPlotEveryLevel, opts.debugPlotPause);

% ===== Mesh plan (used when doDryRun=false) =====
Hmov_list = [4.0 2.0 1.0 0.5 0.25];   % mm (moving/slave mesh sizes coarse->fine)
Hfix_fine = 0.25;                    % mm (fixed/master mesh size)
pushIn    = 0.5 * tol;               % mm (keep <= tol for planar-ish contacts)

% ===== Dry run: placement + plan check without remeshing =====
doDryRun = true;                     % set false once happy with geometry

fprintf('Mesh plan (only used in remesh mode):\n');
fprintf('  Master (fixed) Hmax = %.3f mm\n', Hfix_fine);
fprintf('  Slave  (moving) Hmax levels = %s mm\n', mat2str(Hmov_list));
fprintf('  pushIn = %.4f mm (<= tol=%.4f)\n\n', pushIn, tol);

if doDryRun
    fprintf('=== DRY RUN: placement only (uses original STL triangles, no remeshing) ===\n');

    [Fc, Vc] = loadMeshAny(cubeStl);
    [Fs, Vs] = loadMeshAny(sphereStl);
    [Fg, Vg] = loadMeshAny(grooveStl);

    cubeBody   = buildBodyStruct(Fc, Vc);
    sphereBody = buildBodyStruct(Fs, Vs);
    grooveBody = buildBodyStruct(Fg, Vg);

    fprintf('\n--- Dry run 1/3: cube-cube ---\n');
    [VfixP, VmovP] = placeTwoBodies_overlap(Vc, Vc, Fc, Fc, 'cube', pushIn);
    fixedBody = cubeBody; fixedBody.V = VfixP; fixedBody.bbox = bboxFromV(VfixP);
    movBody   = cubeBody; movBody.V   = VmovP; movBody.bbox   = bboxFromV(VmovP);
    debugPlotPair('cube-cube (dry run)', 1, fixedBody, movBody, tol, opts);
    printPlacementChecks(VfixP, VmovP, 'cube');

    fprintf('\n--- Dry run 2/3: sphere-sphere ---\n');
    [VfixP, VmovP] = placeTwoBodies_overlap(Vs, Vs, Fs, Fs, 'sphere', pushIn);
    fixedBody = sphereBody; fixedBody.V = VfixP; fixedBody.bbox = bboxFromV(VfixP);
    movBody   = sphereBody; movBody.V   = VmovP; movBody.bbox   = bboxFromV(VmovP);
    figure(1002); clf
    debugPlotPair('sphere-sphere (dry run)', 1, fixedBody, movBody, tol, opts);
    printPlacementChecks(VfixP, VmovP, 'sphere');

    fprintf('\n--- Dry run 3/3: sphere-in-groove ---\n');
    [VfixP, VmovP] = placeSphereInGroove(Vg, Vs, Fg, Fs, pushIn);
    fixedBody = grooveBody; fixedBody.V = VfixP; fixedBody.bbox = bboxFromV(VfixP);
    movBody   = sphereBody; movBody.V   = VmovP; movBody.bbox   = bboxFromV(VmovP);
    figure(1003); clf
    debugPlotPair('sphere-in-groove (dry run)', 1, fixedBody, movBody, tol, opts);
    printPlacementChecks(VfixP, VmovP, 'sphere');

    fprintf('\nDry run complete. If placements are correct, set doDryRun=false and rerun.\n');
    return
end

% ===== Full remesh convergence run =====
fprintf('=== REMESH RUN: generateMesh surface remeshing + convergence ===\n');

Results = struct('name',{},'A',{},'Nfaces_moving',{},'dA_pct',{});

fprintf('--- Pair 1/3: cube-cube (refine moving cube only) ---\n');
Results(1) = runPair_generateMesh( ...
    'cube-cube', cubeStl, cubeStl, Hfix_fine, Hmov_list, ...
    @(Vfix,Vmov,Ffix,Fmov) placeTwoBodies_overlap(Vfix,Vmov,Ffix,Fmov,'cube',pushIn), tol, opts);

fprintf('\n--- Pair 2/3: sphere-sphere (refine moving sphere only) ---\n');
Results(2) = runPair_generateMesh( ...
    'sphere-sphere', sphereStl, sphereStl, Hfix_fine, Hmov_list, ...
    @(Vfix,Vmov,Ffix,Fmov) placeTwoBodies_overlap(Vfix,Vmov,Ffix,Fmov,'sphere',pushIn), tol, opts);

fprintf('\n--- Pair 3/3: sphere-in-groove (refine moving sphere only) ---\n');
Results(3) = runPair_generateMesh( ...
    'sphere-in-groove', grooveStl, sphereStl, Hfix_fine, Hmov_list, ...
    @(Vfix,Vmov,Ffix,Fmov) placeSphereInGroove(Vfix,Vmov,Ffix,Fmov,pushIn), tol, opts);

figure; hold on; grid on
xlabel('log10(#faces of moving mesh)')
ylabel('Contact area (mm^2)')
title(sprintf('Mesh convergence (3 pairs), tol = %.3f mm', tol))
for i = 1:numel(Results)
    plot(log10(Results(i).Nfaces_moving), Results(i).A, '-o', 'DisplayName', Results(i).name);
end
legend('Location','best')

figure; hold on; grid on
xlabel('log10(#faces of moving mesh)')
ylabel('% change vs previous level')
title('Convergence rate (% change)')
for i = 1:numel(Results)
    plot(log10(Results(i).Nfaces_moving), Results(i).dA_pct, '-o', 'DisplayName', Results(i).name);
end
legend('Location','best')

fprintf('\n=== Final convergence summary ===\n');
for i = 1:numel(Results)
    r = Results(i);
    fprintf('%-16s | levels=%2d | faces: %s\n', r.name, numel(r.A), mat2str(r.Nfaces_moving));
    fprintf('  areas: %s\n', mat2str(round(r.A,3)));
    fprintf('  dA%% :  %s\n', mat2str(round(r.dA_pct,3)));
    fprintf('  A_final = %.4f mm^2, last %%Δ = %.4f\n\n', r.A(end), r.dA_pct(end));
end

% ============================ FUNCTIONS ============================

function out = runPair_generateMesh(name, fixStl, movStl, Hfix, Hmov_list, placer, tol, opts)
    [Ffix, Vfix] = remeshSurface_generateMesh(fixStl, Hfix);
    fixedBody0 = buildBodyStruct(Ffix, Vfix);
    fixedBody0.kdtree = KDTreeSearcher(fixedBody0.centroids);
    Vfix0 = fixedBody0.V;

    nLevels = numel(Hmov_list);
    A  = nan(1,nLevels);
    Nf = nan(1,nLevels);
    dA = nan(1,nLevels);

    for k = 1:nLevels
        Hm = Hmov_list(k);
        [Fmov, Vmov] = remeshSurface_generateMesh(movStl, Hm);
        Nf(k) = size(Fmov,1);

        fprintf('  Level %d/%d | Hmov=%.3f mm | moving faces=%d ... ', k, nLevels, Hm, Nf(k));

        [VfixP, VmovP] = placer(Vfix0, Vmov, Ffix, Fmov);

        fixedBody = fixedBody0;
        fixedBody.V    = VfixP;
        fixedBody.bbox = bboxFromV(VfixP);

        movBody = buildBodyStruct(Fmov, VmovP);

        zGap = min(VmovP(:,3)) - max(VfixP(:,3));
        fprintf('zGap=%.4f mm | ', zGap);

        if isfield(opts,'debugPlot') && opts.debugPlot
            doPlot = true;
            if isfield(opts,'debugPlotEveryLevel') && ~opts.debugPlotEveryLevel
                doPlot = (k==1);
            end
            if doPlot
                debugPlotPair(name, k, fixedBody, movBody, tol, opts);
                if isfield(opts,'debugPlotPause') && opts.debugPlotPause
                    fprintf('Press any key in figure to continue...\n');
                    waitforbuttonpress;
                end
            end
        end

        res = computeContactArea_STS(fixedBody, movBody, tol, opts);
        if isstruct(res), A(k) = res.contactArea; else, A(k) = res; end

        if k >= 2 && A(k-1) ~= 0
            dA(k) = 100 * abs(A(k) - A(k-1)) / abs(A(k-1));
        end

        if A(k) == 0
            fprintf('A=0 (sanity: increase pushIn a bit, or increase tol)\n');
        else
            if k == 1
                fprintf('A=%.4f\n', A(k));
            else
                fprintf('A=%.4f | dA%%=%.3f\n', A(k), dA(k));
            end
        end
    end

    out = struct();
    out.name = name;
    out.A = A;
    out.Nfaces_moving = Nf;
    out.dA_pct = dA;

    fprintf('  -> %s final: A=%.4f, last dA%%=%.3f\n', name, A(end), dA(end));
end

function [Fsurf, Vsurf] = remeshSurface_generateMesh(stlPath, Hmax)
    m = createpde(1);
    importGeometry(m, stlPath);
    msh = generateMesh(m, 'Hmax', Hmax);

    DT = triangulation(msh.Elements', msh.Nodes');
    [Fsurf, Vsurf] = freeBoundary(DT);

    [Vsurf, ~, ic] = unique(round(Vsurf,12), 'rows', 'stable');
    Fsurf = ic(Fsurf);
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

function [VfixP, VmovP] = placeTwoBodies_overlap(Vfix, Vmov, Ffix, Fmov, shapeType, pushIn)
    VfixP = Vfix;

    cFix = mean(Vfix,1);
    cMov = mean(Vmov,1);
    VmovP = Vmov + (cFix - cMov);

    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(VmovP);

    if strcmpi(shapeType,'sphere')
        Rfix = 0.5 * (bbFix(1,2) - bbFix(1,1));
        Rmov = 0.5 * (bbMov(1,2) - bbMov(1,1));
        targetCenterDist = (Rfix + Rmov) - pushIn;
    else
        Hfix = 0.5 * (bbFix(3,2) - bbFix(3,1));
        Hmov = 0.5 * (bbMov(3,2) - bbMov(3,1));
        targetCenterDist = (Hfix + Hmov) - pushIn;
    end

    cFix = mean(VfixP,1);
    cMov = mean(VmovP,1);

    dz = targetCenterDist - (cMov(3) - cFix(3));
    VmovP(:,3) = VmovP(:,3) + dz;
end

function [VfixP, VmovP] = placeSphereInGroove(Vfix, Vmov, Ffix, Fmov, pushIn)
    VfixP = Vfix;

    cFixXY = mean(Vfix(:,1:2),1);
    cMovXY = mean(Vmov(:,1:2),1);

    VmovP = Vmov;
    VmovP(:,1:2) = VmovP(:,1:2) + (cFixXY - cMovXY);

    bbFix = bboxFromV(VfixP);

    zTargetCenter = 0.5*(bbFix(3,1) + bbFix(3,2));
    cMov = mean(VmovP,1);
    VmovP(:,3) = VmovP(:,3) + (zTargetCenter - cMov(3));

    zFixTop = bbFix(3,2);
    zMovBottom = min(VmovP(:,3));
    dzDown = (zFixTop - pushIn) - zMovBottom;
    VmovP(:,3) = VmovP(:,3) + dzDown;

    VmovP(:,3) = VmovP(:,3) - 0.25*pushIn;
end

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
end

function printPlacementChecks(VfixP, VmovP, type)
    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(VmovP);
    cFix  = mean(VfixP,1);
    cMov  = mean(VmovP,1);

    fprintf('  Fixed bbox Z: [%.3f, %.3f]\n', bbFix(3,1), bbFix(3,2));
    fprintf('  Mov   bbox Z: [%.3f, %.3f]\n', bbMov(3,1), bbMov(3,2));
    fprintf('  Centre distance (XYZ): %.4f mm\n', norm(cMov - cFix));

    if strcmpi(type,'sphere')
        Rfix = 0.5*(bbFix(1,2)-bbFix(1,1));
        Rmov = 0.5*(bbMov(1,2)-bbMov(1,1));
        sep  = norm(cMov - cFix);
        overlap = (Rfix + Rmov) - sep;
        fprintf('  Sphere approx radii: Rfix=%.3f, Rmov=%.3f\n', Rfix, Rmov);
        fprintf('  Approx overlap (Rfix+Rmov - centreDist) = %.4f mm\n', overlap);
    else
        zFixTop = bbFix(3,2);
        zFixBot = bbFix(3,1);
        zMovTop = bbMov(3,2);
        zMovBot = bbMov(3,1);
        zOverlap = min(zFixTop, zMovTop) - max(zFixBot, zMovBot);
        fprintf('  Z-overlap (min tops - max bottoms) = %.4f mm\n', zOverlap);
    end
end

function debugPlotPair(pairName, levelIdx, fixedBody, movBody, tol, opts)
    figure(1001); clf
    set(gcf,'Name',sprintf('%s | level %d', pairName, levelIdx));
    hold on; axis equal; grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s | level %d | tol=%.3f', pairName, levelIdx, tol));

    patch('Faces', fixedBody.F, 'Vertices', fixedBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.35);

    patch('Faces', movBody.F, 'Vertices', movBody.V, ...
        'FaceColor', 'none', 'EdgeColor', [0 0.7 0], 'EdgeAlpha', 0.35);

    if isfield(opts,'plotCentroids') && opts.plotCentroids
        scatter3(fixedBody.centroids(:,1), fixedBody.centroids(:,2), fixedBody.centroids(:,3), 5, '.');
        scatter3(movBody.centroids(:,1), movBody.centroids(:,2), movBody.centroids(:,3), 5, '.');
    end

    roiMin = max(movBody.bbox(:,1), fixedBody.bbox(:,1)) - opts.roiExpandFactor * tol;
    roiMax = min(movBody.bbox(:,2), fixedBody.bbox(:,2)) + opts.roiExpandFactor * tol;
    drawBBox3D(roiMin, roiMax);

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

function varargout = computeContactArea_STS(master, slave, varargin)
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

    if ~isfield(opts,'tol'),               opts.tol               = 0.2; end
    if ~isfield(opts,'roiExpandFactor'),   opts.roiExpandFactor   = 1.5; end
    if ~isfield(opts,'neighRadiusFactor'), opts.neighRadiusFactor = 5.0; end
    if ~isfield(opts,'maxNeighbours'),     opts.maxNeighbours     = 25;  end
    if ~isfield(opts,'sampleMode'),        opts.sampleMode        = 'centroid'; end
    if ~isfield(opts,'sampleThreshold'),   opts.sampleThreshold   = 0.5; end
    if ~isfield(opts,'returnMask'),        opts.returnMask        = (nargout > 1); end
    if ~isfield(opts,'useCentroidPrefilter'), opts.useCentroidPrefilter = true; end
    if ~isfield(opts,'prefilterFactor'),      opts.prefilterFactor      = 6.0;  end

    if ~isempty(tol)
        opts.tol = tol;
    end
    tol = opts.tol;

    master = ensureBodyFields_fast(master);
    slave  = ensureBodyFields_fast(slave);

    nSlaveTris = size(slave.F, 1);

    if ~isfield(master,'kdtree') || isempty(master.kdtree)
        master.kdtree = KDTreeSearcher(master.centroids);
    end

    roiMin = max(slave.bbox(:,1), master.bbox(:,1));
    roiMax = min(slave.bbox(:,2), master.bbox(:,2));
    roiMin = roiMin - opts.roiExpandFactor * tol;
    roiMax = roiMax + opts.roiExpandFactor * tol;

    cS = slave.centroids;
    inROI = cS(:,1) >= roiMin(1) & cS(:,1) <= roiMax(1) & ...
            cS(:,2) >= roiMin(2) & cS(:,2) <= roiMax(2) & ...
            cS(:,3) >= roiMin(3) & cS(:,3) <= roiMax(3);

    candidateIdx = find(inROI);
    if isempty(candidateIdx)
        A_contact = 0;
        if opts.returnMask
            contactMask = false(nSlaveTris,1);
        else
            contactMask = [];
        end
        [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
        return
    end

    neighRadius = max(1e-12, opts.neighRadiusFactor * tol);
    if opts.useCentroidPrefilter
        cCand = cS(candidateIdx,:);
        [~, d0] = knnsearch(master.kdtree, cCand, 'K', 1);
        gate = max(neighRadius, opts.prefilterFactor * tol);
        keep = (d0 <= gate);
        candidateIdx = candidateIdx(keep);

        if isempty(candidateIdx)
            A_contact = 0;
            if opts.returnMask
                contactMask = false(nSlaveTris,1);
            else
                contactMask = [];
            end
            [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
            return
        end
    end

    Fm = master.F; Vm = master.V;
    Am = Vm(Fm(:,1),:);
    Bm = Vm(Fm(:,2),:);
    Cm = Vm(Fm(:,3),:);

    Fs = slave.F;  Vs = slave.V;

    if opts.returnMask
        contactMask = false(nSlaveTris, 1);
    else
        contactMask = [];
    end

    K = max(1, opts.maxNeighbours);
    A_contact = 0;

    for ii = 1:numel(candidateIdx)
        tIdx = candidateIdx(ii);

        sTri = Fs(tIdx,:);
        p1 = Vs(sTri(1),:);
        p2 = Vs(sTri(2),:);
        p3 = Vs(sTri(3),:);

        ps = samplePointsFast(p1,p2,p3, slave.areas(tIdx), opts);

        hit = false;
        for sp = 1:size(ps,1)
            q = ps(sp,:);

            [neighIdx, neighDist] = knnsearch(master.kdtree, q, 'K', K);
            inR = neighDist <= neighRadius;
            if ~any(inR)
                continue
            end
            neighIdx = neighIdx(inR);

            for k = 1:numel(neighIdx)
                j = neighIdx(k);
                d = pointTriangleDistance_fast(q, Am(j,:), Bm(j,:), Cm(j,:));
                if d <= tol
                    hit = true;
                    break
                end
            end
            if hit
                break
            end
        end

        if hit
            A_contact = A_contact + slave.areas(tIdx);
            if opts.returnMask
                contactMask(tIdx) = true;
            end
        end
    end

    [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
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

function ps = samplePointsFast(p1,p2,p3, triArea, opts)
    c  = (p1+p2+p3)/3;
    if strcmpi(opts.sampleMode,'adaptive') && triArea > opts.sampleThreshold
        m12 = (p1+p2)/2;
        m23 = (p2+p3)/2;
        m31 = (p3+p1)/2;
        ps = [c; m12; m23; m31];
    else
        ps = c;
    end
end

function d = pointTriangleDistance_fast(p, a, b, c)
    ab = b - a;
    ac = c - a;
    ap = p - a;

    d1 = dot(ab, ap);
    d2 = dot(ac, ap);
    if d1 <= 0 && d2 <= 0
        d = norm(ap); return;
    end

    bp = p - b;
    d3 = dot(ab, bp);
    d4 = dot(ac, bp);
    if d3 >= 0 && d4 <= d3
        d = norm(bp); return;
    end

    vc = d1*d4 - d3*d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3);
        proj = a + v * ab;
        d = norm(p - proj); return;
    end

    cp = p - c;
    d5 = dot(ab, cp);
    d6 = dot(ac, cp);
    if d6 >= 0 && d5 <= d6
        d = norm(cp); return;
    end

    vb = d5*d2 - d1*d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6);
        proj = a + w * ac;
        d = norm(p - proj); return;
    end

    va = d3*d6 - d5*d4;
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        proj = b + w * (c - b);
        d = norm(p - proj); return;
    end

    denom = 1 / (va + vb + vc);
    v = vb * denom;
    w = vc * denom;
    proj = a + ab * v + ac * w;
    d = norm(p - proj);
end

function varargout = packOutputs(A_contact, contactMask)
    if nargout <= 1
        results = struct();
        results.contactArea = A_contact;
        results.contactMask = contactMask;
        varargout = {results};
    else
        varargout = {A_contact, contactMask};
    end
end
