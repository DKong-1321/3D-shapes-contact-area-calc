% tolerance_sweep_singlePose.m
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir,'model'),'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end
modelDir = fullfile(projectRoot,'model');

planeStlPath  = fullfile(modelDir,'50mm cube 4 iterations.stl');
roundStlPath  = fullfile(modelDir,'50mm diameter sphere.stl');

baseOpts = struct();
baseOpts.sampleThreshold   = 0.5;
baseOpts.neighRadiusFactor = 2.0;
baseOpts.maxNeighbours     = 15;
baseOpts.roiExpandFactor   = 1.05;

errThreshPct = 5;
needConsec   = 3;

tol0        = 0.50;
tolFactor   = 0.80;
tolMin      = 0.02;
maxTolIters = 25;

HmaxChosen  = 0.35;   % set this to whatever you decided from mesh study

[F_plane, V_plane] = loadStlMesh(planeStlPath, 'plane');
master = buildBodyStruct(F_plane, V_plane);

fprintf('\nGenerating sphere mesh once at Hmax = %.3f mm...\n', HmaxChosen);

% Option A (your current route): PDE remesh once
[F_round, V_round] = generateRoundMeshFromSTL(roundStlPath, HmaxChosen);

% Option B (faster if you have a saved mesh):
% load(fullfile(modelDir,'sphere_H0p35.mat'),'F_round','V_round');

Tpose = placeSlaveOnTopOfMasterByAABB(master.V, V_round);

V_round = applyT(V_round, Tpose);
slave   = buildBodyStruct(F_round, V_round);

fprintf('Placement translation: [%.3f, %.3f, %.3f] mm\n', Tpose(1,4), Tpose(2,4), Tpose(3,4));
fprintf('Sphere faces: %d\n\n', size(F_round,1));

tolRows = [];
tol = tol0;
consecGood = 0;
A_prev = NaN;

fprintf('=== Tolerance sweep (single pose) ===\n');

for iter = 1:maxTolIters
    if tol < tolMin
        fprintf('Stopping: tol < tolMin (%.4f < %.4f)\n', tol, tolMin);
        break;
    end

    fprintf('Tol iter %d/%d | tol = %.4f mm\n', iter, maxTolIters, tol);

    opts = baseOpts;
    opts.tol = tol;

    % cheap skip: if clearly separated at this tol, area is zero
    if aabbGap(master.V, slave.V) > opts.tol
        A = 0;
    else
        results = computeContactArea_STS(master, slave, opts);
        A = results.contactArea;
    end

    if isnan(A_prev)
        dA_pct = NaN;
    else
        dA_pct = 100 * abs(A - A_prev) / max(eps, abs(A_prev));
    end

    tolRows = [tolRows; 1, iter, tol, A, dA_pct]; %#ok<AGROW>

    if ~isnan(dA_pct) && dA_pct < errThreshPct
        consecGood = consecGood + 1;
    else
        consecGood = 0;
    end

    if consecGood >= needConsec
        fprintf('  ✓ Converged after %d tolerance steps\n', iter);
        break;
    end

    A_prev = A;
    tol = tol * tolFactor;
end

tolT = array2table(tolRows, 'VariableNames', ...
    {'PosID','TolIter','Tol_mm','ContactArea_mm2','DeltaArea_pct'});
disp(tolT);

figure;
plot(tolT.Tol_mm, tolT.ContactArea_mm2, '-o', 'LineWidth', 1.5);
set(gca,'XDir','reverse');
xlabel('Tolerance (mm)'); ylabel('Contact area (mm^2)'); grid on;

fprintf('\nDone.\n');

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

function [F_surf, V_surf] = generateRoundMeshFromSTL(stlPath, Hmax)
    model = femodel(Geometry=stlPath);
    model = generateMesh(model, "Hmax", Hmax, "GeometricOrder", "linear");
    msh = model.Geometry.Mesh;
    TR  = triangulation(msh.Elements', msh.Nodes');
    [F_surf, V_surf] = freeBoundary(TR);
    [F_surf, V_surf] = patchCleanUnused(F_surf, V_surf);
end

function Vt = applyT(V, T)
    Vh = [V, ones(size(V,1),1)];
    Vt = (T * Vh')';
    Vt = Vt(:,1:3);
end

function [F2,V2] = patchCleanUnused(F,V)
    used = unique(F(:));
    map  = zeros(max(used),1);
    map(used) = 1:numel(used);
    V2 = V(used,:);
    F2 = map(F);
end

function g = aabbGap(Va, Vb)
    aMin = min(Va,[],1); aMax = max(Va,[],1);
    bMin = min(Vb,[],1); bMax = max(Vb,[],1);
    dx = max([bMin(1)-aMax(1), aMin(1)-bMax(1), 0]);
    dy = max([bMin(2)-aMax(2), aMin(2)-bMax(2), 0]);
    dz = max([bMin(3)-aMax(3), aMin(3)-bMax(3), 0]);
    g = sqrt(dx*dx + dy*dy + dz*dz);
end
