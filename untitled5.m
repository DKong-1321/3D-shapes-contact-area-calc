% sphereSphere_toleranceEdgeOnly_FINE_FINE_demo.m
% -------------------------------------------------------------------------
% DEMO (what you asked):
% - BOTH spheres are the SAME *refined* STL (no generateMesh / no coarse mesh)
% - Place them with a specified penetration (mm)
% - Run the TOLERANCE-ONLY contact method (computeContactArea_STS)
% - Plot shows that under penetration, tolerance contact tends to mark
%   mainly a thin band around the edge of the overlap region (limitation demo)
%
% REQUIREMENTS on path (src/):
%   - buildBodyStruct
%   - computeContactArea_STS   (your tolerance function)
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

fprintf('\n=== Sphere–sphere tolerance contact | FINE vs FINE | edge-only under penetration demo ===\n');
fprintf('Project root: %s\n', projectRoot);
fprintf('Model dir:    %s\n\n', modelDir);

%% ---- FILE (EDIT TO YOUR EXACT FILENAME) ----
sphereStlFine = fullfile(modelDir,'Sphere 50mm diameter 3 it.stl');
assert(exist(sphereStlFine,'file')==2, 'Sphere STL not found: %s', sphereStlFine);

%% ---- SETTINGS ----
penetration_mm = 5.0;   % 0 = just-touching, >0 = overlap/penetration
tol_mm         = 0.20;  % tolerance for computeContactArea_STS
% If the magenta band is too thin to see, try tol_mm = 0.4 or 0.6 (for demo only)

showEdgesFixed  = false;
showEdgesMoving = true;

%% ---- LOAD BOTH SPHERES (FINE) ----
[Ffix, Vfix] = loadMeshAny(sphereStlFine);
[Fmov, Vmov] = loadMeshAny(sphereStlFine);

fixed0  = buildBodyStruct(Ffix, Vfix);
moving0 = buildBodyStruct(Fmov, Vmov);

%% ---- PLACE SPHERE-SPHERE (STACKED +Z, then penetrated) ----
[VfixP, VmovP] = placeTwoSpheres_stackZ(fixed0.V, moving0.V, penetration_mm);

fixed  = buildBodyStruct(fixed0.F,  VfixP);
moving = buildBodyStruct(moving0.F, VmovP);

%% ---- RUN TOLERANCE CONTACT ----
opts = struct();
opts.tol = tol_mm;

% (match your computeContactArea_STS defaults, but explicit = reproducible)
opts.roiExpandFactor      = 1.5;
opts.neighRadiusFactor    = 5.0;
opts.maxNeighbours        = 25;
opts.sampleMode           = 'centroid';  % keep this to demonstrate the limitation clearly
opts.sampleThreshold      = 0.5;
opts.returnMask           = true;
opts.useCentroidPrefilter = true;
opts.prefilterFactor      = 6.0;

res  = computeContactArea_STS(fixed, moving, tol_mm, opts);
mask = res.contactMask;
A    = res.contactArea;

fprintf('Requested penetration: %.3f mm\n', penetration_mm);
fprintf('Tolerance:             %.3f mm\n', tol_mm);
fprintf('Tolerance-contact area: %.6f mm^2\n', A);
fprintf('Marked triangles:       %d / %d\n\n', nnz(mask), size(moving.F,1));

%% ---- DIAGNOSTIC: measured penetration from bbox radii (sanity) ----
bbFix = bboxFromV(VfixP);
bbMov = bboxFromV(VmovP);
Rfix = 0.5*(bbFix(1,2)-bbFix(1,1));
Rmov = 0.5*(bbMov(1,2)-bbMov(1,1));
cFix = bboxCenter(VfixP);
cMov = bboxCenter(VmovP);
pen_meas = (Rfix + Rmov) - norm(cMov - cFix);
fprintf('Measured penetration (bbox-based): %.6f mm\n\n', pen_meas);

%% ---- PLOT: EXPECT "EDGE BAND" UNDER PENETRATION ----
figure('Color','w'); hold on; axis equal; grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Tolerance-only contact (fine/fine) | pen=%.2f mm | tol=%.2f mm | A=%.3f mm^2', ...
    penetration_mm, tol_mm, A), 'Interpreter','none');

% Fixed sphere (grey)
patch('Faces', fixed.F, 'Vertices', fixed.V, ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.08, ...
    'EdgeColor', tern(showEdgesFixed,[0 0 0],'none'), ...
    'EdgeAlpha', 0.12, 'LineWidth', 0.20);

% Moving sphere (green)
patch('Faces', moving.F, 'Vertices', moving.V, ...
    'FaceColor', [0 0.7 0], 'FaceAlpha', 0.08, ...
    'EdgeColor', tern(showEdgesMoving,[0 0 0],'none'), ...
    'EdgeAlpha', 0.15, 'LineWidth', 0.20);

% Tolerance contact triangles (magenta) — should appear as a thin ring/band
if any(mask)
    patch('Faces', moving.F(mask,:), 'Vertices', moving.V, ...
        'FaceColor', [1 0 1], 'FaceAlpha', 0.95, ...
        'EdgeColor', [0 0 0], 'EdgeAlpha', 0.20, ...
        'LineWidth', 0.20);
else
    text(0,0,0,'No triangles marked by tolerance method (try larger tol).','FontSize',12);
end

view(3); camlight; lighting gouraud

%% ============================ HELPERS ============================

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
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

function [VfixP, VmovP] = placeTwoSpheres_stackZ(Vfix, Vmov, penetration)
% Align XY by bbox centres, stack moving ABOVE fixed along +Z,
% then reduce centre distance by "penetration" to create overlap.

    VfixP = Vfix;

    cFix = bboxCenter(VfixP);
    cMov = bboxCenter(Vmov);

    bbFix = bboxFromV(VfixP);
    bbMov = bboxFromV(Vmov);
    Rfix = 0.5 * (bbFix(1,2) - bbFix(1,1));
    Rmov = 0.5 * (bbMov(1,2) - bbMov(1,1));

    penetration = max(0, min(penetration, (Rfix + Rmov)));

    % Align XY
    VmovP = Vmov;
    VmovP(:,1) = VmovP(:,1) + (cFix(1) - cMov(1));
    VmovP(:,2) = VmovP(:,2) + (cFix(2) - cMov(2));

    % Recompute centre after XY shift
    cMov2 = bboxCenter(VmovP);

    targetCenterDist = (Rfix + Rmov) - penetration;

    % Stack above (+Z)
    dz = (cFix(3) + targetCenterDist) - cMov2(3);
    VmovP(:,3) = VmovP(:,3) + dz;
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
