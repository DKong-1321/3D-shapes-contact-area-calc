% toleranceTable_autoStop_fullSphere_showContact_RED.m
% Builds a tolerance table by *reducing tolerance step-by-step* until the
% % change in contact area between consecutive tolerances is < errThreshPct
% for needConsec consecutive steps (then stops).
% Also (optionally) saves a separate figure for each tolerance showing the
% FULL sphere + highlighted contact faces (RED), using ROI slab for speed.



%% ===== PATHS =====
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
modelDir = fullfile(projectRoot,'model');

cubeStlPath   = fullfile(modelDir,'50mm cube 4 iterations.stl');
sphereStlPath = fullfile(modelDir,'50mm diameter sphere.stl');  % already refined STL

figDir = scriptDir;
framesDir = fullfile(figDir,'tol_frames');
if ~exist(framesDir,'dir'), mkdir(framesDir); end

%% ===== SETTINGS =====
baseOpts = struct();
baseOpts.sampleThreshold   = 0.5;
baseOpts.neighRadiusFactor = 2.0;
baseOpts.maxNeighbours     = 15;
baseOpts.roiExpandFactor   = 1.05;

% --- tolerance schedule (start high, reduce each step) ---
tol0        = 0.50;
tolFactor   = 0.80;
tolMin      = 0.02;
maxTolIters = 25;

% --- auto-stop criterion (this is the "tolerance table idea") ---
errThreshPct = 5;     % stop when ΔA% < this...
needConsec   = 3;     % ...for this many consecutive tolerance reductions

useSlabROI  = true;   % recommended for speed

saveEveryStepPng   = true;  % per tolerance step figure saved
showFiguresOnScreen = false; % set true if you want tabs popping up

% Print fields once to diagnose whether contact face indices exist
printResultsFieldnamesOnce = false;

%% ===== LOAD =====
[F_cube, V_cube] = loadStlMesh(cubeStlPath, 'plane');
master = buildBodyStruct(F_cube, V_cube);
zTop = max(master.V(:,3));

[F_sph0, V_sph0] = loadStlMesh(sphereStlPath, 'round');
fprintf('Loaded sphere STL faces: %d\n', size(F_sph0,1));

%% ===== PLACE SPHERE ON TOP OF CUBE (MATRIX) =====
Tpose = placeSlaveOnTopOfMasterByAABB(master.V, V_sph0);
V_sph = applyT(V_sph0, Tpose);
fprintf('Placement translation: [%.3f, %.3f, %.3f] mm\n', Tpose(1,4), Tpose(2,4), Tpose(3,4));

fullSlave = buildBodyStruct(F_sph0, V_sph); % full sphere body for display

%% ===== SWEEP UNTIL CONVERGED =====
tolVals         = nan(maxTolIters,1);
Acontact        = nan(maxTolIters,1);
dA_pct          = nan(maxTolIters,1);
ROIFaces        = nan(maxTolIters,1);
NumContactFaces = nan(maxTolIters,1);

tol = tol0;
k = 0;
consecGood = 0;

fprintf('\n=== Auto-stop tolerance sweep (FULL sphere shown + RED contact) ===\n');
fprintf('Stop rule: ΔA < %.1f%% for %d consecutive tol reductions\n', errThreshPct, needConsec);

while true
    if k >= maxTolIters
        warning('Reached maxTolIters without satisfying convergence criterion.');
        break;
    end
    if tol < tolMin
        warning('Reached tolMin without satisfying convergence criterion.');
        break;
    end

    k = k + 1;
    fprintf('Tol iter %d/%d | tol = %.4f mm\n', k, maxTolIters, tol);

    opts = baseOpts;
    opts.tol = tol;

    contactIdxFull = [];
    A = 0;
    nROI = 0;

    if useSlabROI
        triIdxROI_full = roiByTopSlab(V_sph, F_sph0, zTop, tol, baseOpts.roiExpandFactor); % FULL face indices
        nROI = numel(triIdxROI_full);

        if ~isempty(triIdxROI_full)
            [F_roi, V_roi, backMap] = extractSubPatchWithBackMap(F_sph0, V_sph, triIdxROI_full);
            slaveROI = buildBodyStruct(F_roi, V_roi);

            results = computeContactArea_STS(master, slaveROI, opts);
            A = results.contactArea;

            if printResultsFieldnamesOnce
                fprintf('--- fieldnames(results) ---\n');
                disp(fieldnames(results));
                if isfield(results,'contact') && isstruct(results.contact)
                    fprintf('--- fieldnames(results.contact) ---\n');
                    disp(fieldnames(results.contact));
                end
                printResultsFieldnamesOnce = false;
            end

            contactIdxROI  = getSlaveContactFaceIdx(results);      % ROI face indices
            contactIdxFull = mapIdxToFull(contactIdxROI, backMap); % mapped to FULL indices
        end
    else
        results = computeContactArea_STS(master, fullSlave, opts);
        A = results.contactArea;

        if printResultsFieldnamesOnce
            fprintf('--- fieldnames(results) ---\n');
            disp(fieldnames(results));
            if isfield(results,'contact') && isstruct(results.contact)
                fprintf('--- fieldnames(results.contact) ---\n');
                disp(fieldnames(results.contact));
            end
            printResultsFieldnamesOnce = false;
        end

        contactIdxFull = getSlaveContactFaceIdx(results);
        nROI = size(F_sph0,1);
    end

    fprintf('  ROI faces: %d | Contact faces detected: %d | A=%.4f\n', nROI, numel(contactIdxFull), A);

    tolVals(k)         = tol;
    Acontact(k)        = A;
    ROIFaces(k)        = nROI;
    NumContactFaces(k) = numel(contactIdxFull);

    if k == 1
        dA_pct(k) = NaN;
        consecGood = 0;
    else
        dA = abs(Acontact(k) - Acontact(k-1));
        dA_pct(k) = 100 * dA / max(eps, abs(Acontact(k-1)));

        if dA_pct(k) < errThreshPct
            consecGood = consecGood + 1;
        else
            consecGood = 0;
        end
    end

    % Per-step figure (FULL sphere + RED contact)
    if saveEveryStepPng || showFiguresOnScreen
        if showFiguresOnScreen
            fig = figure(100 + k); clf(fig);
            set(fig,'Color','w');
        else
            fig = figure('Visible','off');
            set(fig,'Color','w');
        end

        ttl = sprintf('Step %d | tol=%.4f mm | A=%.3f mm^2 | ΔA=%.2f%% | ROI=%d | contact=%d', ...
            k, tol, A, dA_pct(k), nROI, numel(contactIdxFull));

        renderFullSphereWithContact(fig, master, fullSlave, contactIdxFull, ttl);

        if saveEveryStepPng
            outPng = fullfile(framesDir, sprintf('tolStep_%02d_tol_%0.4f.png', k, tol));
            exportgraphics(fig, outPng, 'Resolution', 200);
        end

        if ~showFiguresOnScreen
            close(fig);
        end
    end

    % Convergence stop
    if consecGood >= needConsec
        fprintf('\nConverged: ΔA < %.1f%% for %d consecutive tol reductions.\n', errThreshPct, needConsec);
        break;
    end

    tol = tol * tolFactor;
end

%% ===== TRIM + TABLE =====
tolVals         = tolVals(1:k);
Acontact        = Acontact(1:k);
dA_pct          = dA_pct(1:k);
ROIFaces        = ROIFaces(1:k);
NumContactFaces = NumContactFaces(1:k);

T = table((1:k)', tolVals, Acontact, dA_pct, ROIFaces, NumContactFaces, ...
    'VariableNames', {'Iter','Tol_mm','ContactArea_mm2','DeltaArea_pct','ROIFaces','NumContactFaces'});
disp(T);

% Recommended tolerance = the tolerance at which the criterion is first satisfied
% i.e. the last tolerance in the run (since we stop immediately when consecGood hits needConsec)
tolRecommended = tolVals(end);
fprintf('Recommended tol (auto-stop): %.4f mm\n', tolRecommended);

%% ===== SUMMARY PLOTS =====
figure; set(gcf,'Color','w');
plot(T.Tol_mm, T.ContactArea_mm2, '-o', 'LineWidth', 1.5);
set(gca,'XDir','reverse');
xlabel('Tolerance (mm)'); ylabel('Contact area (mm^2)'); grid on;
exportgraphics(gcf, fullfile(figDir,'toleranceTable_contactArea_vs_tol.png'), 'Resolution', 300);

figure; set(gcf,'Color','w');
plot(T.Tol_mm, T.DeltaArea_pct, '-o', 'LineWidth', 1.5);
yline(errThreshPct, '--');
set(gca,'XDir','reverse');
xlabel('Tolerance (mm)'); ylabel('\Delta Contact area (%)'); grid on;
exportgraphics(gcf, fullfile(figDir,'toleranceTable_deltaArea_vs_tol.png'), 'Resolution', 300);

fprintf('\nSaved per-step frames in:\n  %s\n', framesDir);
fprintf('Saved plots in:\n  %s\n  %s\n', ...
    fullfile(figDir,'toleranceTable_contactArea_vs_tol.png'), ...
    fullfile(figDir,'toleranceTable_deltaArea_vs_tol.png'));
fprintf('Done.\n');

%% ===== HELPERS =====

function triIdx = roiByTopSlab(V, F, zTop, tol, roiExpandFactor)
    h = tol * roiExpandFactor;
    zTri = reshape(V(F(:),3), size(F,1), 3);
    zMin = min(zTri,[],2);
    zMax = max(zTri,[],2);
    triIdx = find(zMin <= (zTop + h) & zMax >= (zTop - h));
end

function [Fsub, Vsub, backMap] = extractSubPatchWithBackMap(F, V, triIdxFull)
    backMap = triIdxFull(:);   % ROI face i corresponds to FULL face backMap(i)
    FsubFull = F(triIdxFull,:);

    used = unique(FsubFull(:));
    mapV = zeros(max(used),1);
    mapV(used) = 1:numel(used);

    Vsub = V(used,:);
    Fsub = mapV(FsubFull);
end

function idxFull = mapIdxToFull(idxROI, backMap)
    idxFull = [];
    if isempty(idxROI) || isempty(backMap), return; end
    idxROI = idxROI(:);
    idxROI = idxROI(idxROI >= 1 & idxROI <= numel(backMap));
    idxFull = backMap(idxROI);
end

function T = placeSlaveOnTopOfMasterByAABB(Vmaster, Vslave)
    mMin = min(Vmaster,[],1); mMax = max(Vmaster,[],1);
    sMin = min(Vslave,[],1);  sMax = max(Vslave,[],1);
    mCtr = 0.5*(mMin + mMax);
    sCtr = 0.5*(sMin + sMax);
    dx = mCtr(1) - sCtr(1);
    dy = mCtr(2) - sCtr(2);
    dz = mMax(3) - sMin(3);
    T = eye(4);
    T(1,4) = dx; T(2,4) = dy; T(3,4) = dz;
end

function Vt = applyT(V, T)
    Vh = [V, ones(size(V,1),1)];
    Vt = (T * Vh')';
    Vt = Vt(:,1:3);
end

function idx = getSlaveContactFaceIdx(results)
    idx = [];

    candTop = {'slaveContactFaces','contactFacesSlave','contactFaceIdxSlave', ...
               'slaveFaceIdx','slaveTriIdx','contactIdxSlave','contactFaces'};
    for k = 1:numel(candTop)
        f = candTop{k};
        if isfield(results, f) && ~isempty(results.(f))
            idx = results.(f);
            return;
        end
    end

    if isfield(results,'contact') && isstruct(results.contact)
        candNested = {'slaveFaceIdx','slaveTriIdx','facesSlave','contactFacesSlave','idxSlave'};
        for k = 1:numel(candNested)
            f = candNested{k};
            if isfield(results.contact, f) && ~isempty(results.contact.(f))
                idx = results.contact.(f);
                return;
            end
        end
    end
end

function renderFullSphereWithContact(fig, master, fullSlave, contactIdxFull, ttl)
    figure(fig); clf(fig); hold on;

    patch('Faces', master.F, 'Vertices', master.V, ...
        'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.15);

    patch('Faces', fullSlave.F, 'Vertices', fullSlave.V, ...
        'FaceColor', [0.75 0.75 0.75], ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'LineWidth', 0.25, ...
        'FaceAlpha', 0.60);

    if ~isempty(contactIdxFull)
        contactIdxFull = contactIdxFull(:);
        contactIdxFull = contactIdxFull(contactIdxFull >= 1 & contactIdxFull <= size(fullSlave.F,1));
        if ~isempty(contactIdxFull)
            patch('Faces', fullSlave.F(contactIdxFull,:), ...
                  'Vertices', fullSlave.V, ...
                  'FaceColor', [1.0 0.0 0.0], ...
                  'EdgeColor', [0 0 0], ...
                  'LineWidth', 0.8, ...
                  'FaceAlpha', 1.0);
        end
    end

    axis equal; grid on;
    view(3);
    camlight headlight;
    lighting flat;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(ttl, 'Interpreter', 'none');
    drawnow;
end
