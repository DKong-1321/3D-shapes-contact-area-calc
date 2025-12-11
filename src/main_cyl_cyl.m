% main_cyl_on_cyl_tolSweep.m
% Static cylinder-on-cylinder contact; sweep tolerance from 0.05 to 0.5.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end

addpath(genpath(fullfile(projectRoot, 'src')));

modelDir = fullfile(projectRoot, 'model');
cylFile  = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');

% Contact options
opts.sampleThreshold   = 0.5;
opts.neighRadiusFactor = 2.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

fprintf('Loading cylinder STL twice (bottom + top)...\n');
[Fbot, Vbot] = loadStlMesh(cylFile, 'bottom');
[Ftop, Vtop] = loadStlMesh(cylFile, 'top');

fprintf('Mesh loaded: %d vertices, %d faces\n', size(Vbot,1), size(Fbot,1));

% Fixed stacked geometry (no sweep in translation, only tol)
[Vtop_aligned, shift] = placeCylinderOnCylinder(Vtop, Vbot);
fprintf('Top cylinder translated by [%.4f  %.4f  %.4f]\n', shift);

bottom = buildBodyStruct(Fbot, Vbot);
top    = buildBodyStruct(Ftop, Vtop_aligned);

fprintf('Bottom total area : %.3f\n', sum(bottom.triArea));
fprintf('Top total area    : %.3f\n', sum(top.triArea));

master = bottom;
slave  = top;

fprintf('Master surface: bottom cylinder\n');
fprintf('Slave surface : top cylinder\n');

% Build KD-tree once on master (tri centroids as in your original code)
fprintf('Building KD-tree on master centroids...\n');
master.kdtree = KDTreeSearcher(master.triCentroid);

% Tolerance sweep
tolList = 0.05:0.05:0.5;
nTol    = numel(tolList);
A_list  = zeros(nTol,1);
frac_list = zeros(nTol,1);
mask_store = false(size(slave.F,1), nTol);  % store masks if you want to inspect later

fprintf('Sweeping tolerance from %.2f to %.2f...\n', tolList(1), tolList(end));

for k = 1:nTol
    tol = tolList(k);
    fprintf('  tol = %.3f ... ', tol);
    [A_k, mask_k] = computeContactArea_STS(slave, master, tol, opts);

    A_list(k)     = A_k;
    frac_list(k)  = nnz(mask_k) / numel(mask_k);
    mask_store(:,k) = mask_k;

    fprintf('A = %.4f, frac = %.4f\n', A_k, frac_list(k));
end

% Plot contact area vs tolerance
figure('Color','w'); 
plot(tolList, A_list, 'o-','LineWidth',1.5);
xlabel('Tolerance'); ylabel('Contact area');
title('Cylinder-on-cylinder: contact area vs tolerance');
grid on;

% Visualise contact for one tolerance (e.g. the largest tol)
[~, idxShow] = max(tolList);  % or pick a specific index
tolShow  = tolList(idxShow);
maskShow = mask_store(:, idxShow);

fprintf('\nVisualising contact for tol = %.3f\n', tolShow);

figure('Color','w'); hold on; axis equal;
title(sprintf('Cylinder-on-cylinder (tol = %.3f)', tolShow), 'Interpreter','none');
xlabel X; ylabel Y; zlabel Z;

% Master in light grey
patch('Faces',master.F,'Vertices',master.V, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'EdgeColor','none', ...
      'FaceAlpha',0.3);

% Slave: blue by default, red where in contact
contactColors = repmat([0.0 0.8 1.0], size(slave.F,1), 1);      % blue/cyan
contactColors(maskShow,:) = repmat([1.0 0.1 0.1], nnz(maskShow), 1); % red

patch('Faces',slave.F,'Vertices',slave.V, ...
      'FaceVertexCData',contactColors, ...
      'FaceColor','flat', ...
      'EdgeColor','k', ...
      'FaceAlpha',0.9);

legend({'Bottom (master)','Top (slave: red = contact)'});
view(3); camlight; lighting gouraud; grid on;

if exist('showContactOnly','file') == 2
    showContactOnly(slave, maskShow, tolShow);
end
