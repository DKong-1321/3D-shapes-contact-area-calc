% toleranceSweep_penetrationArea.m
% -------------------------------------------------------------------------
% Tolerance sensitivity sweep for the CLIPPED penetration-area method
% using a *fixed, converged moving mesh* (no further refinement).
%
% INPUT EXPECTATION:
%   This script loads a .mat file produced by your mesh convergence run that
%   contains a variable named `Results` (struct array), where each element has:
%     - name
%     - F_master, V_masterPlaced
%     - F_final,  V_finalPlaced
%
% OUTPUT:
%   - Plots: A(tol) and relative error vs A_true
%   - CSV + MAT saved in the current folder
%
% IMPORTANT:
%   This uses the same contact area function you embedded:
%     computeContactArea_STS_penetration (included at bottom here)
% -------------------------------------------------------------------------



%% ===================== USER SETTINGS =====================
% 1) Load converged meshes from your mesh convergence script output
%    (Change the filename to whatever you saved.)
resultsMatFile = 'meshConvergenceResults.mat';  % <-- EDIT (must contain `Results`)

% 2) Tolerance sweep list (mm)
tolList = [0.01 0.02 0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.75 1.00];

% 3) Define what "true" area means:
%    'minTol'  -> A_true = A(tolList smallest)
%    'medianK' -> A_true = median of the first K smallest tolerances
trueMode = 'medianK';
Ktrue = 3;  % only used if trueMode='medianK'

% 4) Optional: also return contact mask per tol (can be big)
storeMasks = false;

% 5) Plot settings
makePlots = true;

%% ===================== LOAD RESULTS =====================
assert(exist(resultsMatFile,'file')==2, 'Cannot find: %s', resultsMatFile);
S = load(resultsMatFile);
assert(isfield(S,'Results'), 'MAT file must contain variable `Results`.');
Results = S.Results;

assert(isstruct(Results) && ~isempty(Results), 'Results is empty/invalid.');

fprintf('\n=== Tolerance sweep loaded: %s ===\n', resultsMatFile);
fprintf('Pairs in Results: %d\n', numel(Results));
for i=1:numel(Results)
    fprintf('  %d) %s\n', i, Results(i).name);
end

%% ===================== RUN SWEEP =====================
opts = struct(); % passed to contact function (not used heavily, but consistent)

Sweep = struct([]);
row = 0;

for p = 1:numel(Results)
    R = Results(p);

    % Build master/slave bodies for this pair using converged meshes (FIXED)
    masterBody = buildBodyStruct(R.F_master, R.V_masterPlaced);
    masterBody.bbox = bboxFromV(masterBody.V);
    masterBody.faceNormals = faceNormalsFromFV(masterBody.F, masterBody.V);
    masterBody.kdtreeFaces = KDTreeSearcher(masterBody.centroids);

    slaveBody  = buildBodyStruct(R.F_final,  R.V_finalPlaced);
    slaveBody.bbox = bboxFromV(slaveBody.V);

    A = nan(size(tolList));
    if storeMasks
        masks = cell(size(tolList));
    else
        masks = [];
    end

    for t = 1:numel(tolList)
        tol = tolList(t);
        opts.tol = tol;

        res = computeContactArea_STS_penetration(masterBody, slaveBody, tol, opts);
        A(t) = res.contactArea;

        if storeMasks
            masks{t} = res.contactMask;
        end

        fprintf('%-16s | tol=%8.4f mm | A=%.6f mm^2\n', R.name, tol, A(t));
    end

    % Decide A_true for this pair
    [tolSorted, ix] = sort(tolList, 'ascend');
    Asorted = A(ix);

    switch lower(trueMode)
        case 'mintol'
            A_true = Asorted(1);
        case 'mediank'
            k = min(Ktrue, numel(Asorted));
            A_true = median(Asorted(1:k));
        otherwise
            error('Unknown trueMode: %s', trueMode);
    end

    relErrPct = nan(size(A));
    if abs(A_true) > 1e-12
        relErrPct = 100 * abs(A - A_true) / abs(A_true);
    end

    % Store in a flat table-like struct array for easy CSV
    for t = 1:numel(tolList)
        row = row + 1;
        Sweep(row).pair = string(R.name);
        Sweep(row).tol_mm = tolList(t);
        Sweep(row).A_mm2  = A(t);
        Sweep(row).A_true_mm2 = A_true;
        Sweep(row).relErr_pct = relErrPct(t);
    end

    % Optional plots per pair
    if makePlots
        figure('Color','w'); grid on; hold on;
        plot(tolList, A, '-o');
        xlabel('tol (mm)');
        ylabel('Penetration contact area A (mm^2)');
        title(sprintf('%s | A(tol) (fixed converged mesh)', R.name), 'Interpreter','none');

        yline(A_true, '--', sprintf('A_{true}=%.6f', A_true), 'LabelVerticalAlignment','bottom');

        figure('Color','w'); grid on; hold on;
        plot(tolList, relErrPct, '-o');
        xlabel('tol (mm)');
        ylabel('Relative error vs A_{true} (%)');
        title(sprintf('%s | tolerance sensitivity', R.name), 'Interpreter','none');
    end

    % Save per-pair detailed mat (optional)
    pairOut(p).name = R.name; %#ok<SAGROW>
    pairOut(p).tolList = tolList; %#ok<SAGROW>
    pairOut(p).A = A; %#ok<SAGROW>
    pairOut(p).A_true = A_true; %#ok<SAGROW>
    pairOut(p).relErrPct = relErrPct; %#ok<SAGROW>
    if storeMasks
        pairOut(p).masks = masks; %#ok<SAGROW>
    end
end

%% ===================== SAVE CSV + MAT =====================
T = struct2table(Sweep);

timestamp = datestr(now,'yyyymmdd_HHMMSS');
csvName = sprintf('toleranceSweep_%s.csv', timestamp);
matName = sprintf('toleranceSweep_%s.mat', timestamp);

writetable(T, csvName);
save(matName, 'Results', 'tolList', 'trueMode', 'Ktrue', 'pairOut', 'T');

fprintf('\n=== Saved ===\n  %s\n  %s\n', csvName, matName);

%% =====================================================================
%%                         HELPER FUNCTIONS
%% =====================================================================

function bb = bboxFromV(V)
    bb = [min(V(:,1)) max(V(:,1)); ...
          min(V(:,2)) max(V(:,2)); ...
          min(V(:,3)) max(V(:,3))];
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
    results.contactArea_surfaceOnly = 0;
    results.contactArea_penetrationOnly = Apen;
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
            poly = [poly; Pk]; %#ok<AGROW>
        end

        if (dk <= 0 && dk2 > 0) || (dk > 0 && dk2 <= 0)
            t = dk / (dk - dk2);
            Pi = Pk + t*(Pk2 - Pk);
            poly = [poly; Pi]; %#ok<AGROW>
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
