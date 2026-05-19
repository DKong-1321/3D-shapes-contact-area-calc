% main_cube_cylinder_singlePos_singleTol.m
% Cube (MASTER) – cylinder (SLAVE) contact at ONE position and ONE tolerance

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end

addpath(genpath(fullfile(projectRoot, 'src')));
modelDir = fullfile(projectRoot, 'model');

cubeFile = fullfile(modelDir, '50mm cube 4 iterations.stl');
cylFile  = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');

% ---- Single tolerance ----
tol = 0.10; % mm

% ---- Contact options ----
opts.roiExpandFactor        = 1.5;
opts.neighRadiusFactor      = 5.0;
opts.maxNeighbours          = 30;
opts.sampleMode             = 'centroid';   % keep fast
opts.sampleThreshold        = 0.2;          % only used if adaptive
opts.useCentroidPrefilter   = true;
opts.prefilterFactor        = 6.0;
opts.returnMask             = true;         % IMPORTANT for colouring

% ---- Load meshes ----
[Fc, Vc] = loadStlMesh(cubeFile, 'cube');
[Fs, Vs] = loadStlMesh(cylFile,  'cylinder');

% ---- Place cylinder on cube (your helper defines the ONE pose) ----
[Vs_aligned, ~] = placeCylinderOnCube(Vs, Vc);

% ---- Build body structs ----
master = buildBodyStruct(Fc, Vc);         % cube MASTER
slave  = buildBodyStruct(Fs, Vs_aligned); % cylinder SLAVE

fprintf('Cube (master):     %d faces, %d verts\n', size(Fc,1), size(Vc,1));
fprintf('Cylinder (slave):  %d faces, %d verts\n', size(Fs,1), size(Vs_aligned,1));

cubeEdge = computeAverageEdgeLength(Fc, Vc);
cylEdge  = computeAverageEdgeLength(Fs, Vs_aligned);
fprintf('Average edge length (cube):     %.4f mm\n', cubeEdge);
fprintf('Average edge length (cylinder): %.4f mm\n', cylEdge);

% ---- Build KD-tree ONCE on master centroids (as your function expects) ----
% buildBodyStruct might store triCentroid; your function maps it -> centroids if needed,
% but the KD-tree is built from master.centroids, so we ensure it's present.
if isfield(master,'triCentroid') && ~isfield(master,'centroids')
    master.centroids = master.triCentroid;
end
master.kdtree = KDTreeSearcher(master.centroids);

% ---- Compute contact (IMPORTANT: master first, slave second) ----
[A_contact, contactMask] = computeContactArea_STS(master, slave, tol, opts);

fprintf('\nTol = %.3f mm\n', tol);
fprintf('Contact area = %.6f mm^2\n', A_contact);
fprintf('Contact tris (on cylinder) = %d / %d\n', nnz(contactMask), size(Fs,1));

% ---- Plot ----
cubeColor    = [0.8 0.8 0.8];
cylBaseColor = [0.92 0.92 0.92];
contactColor = [0.0 1.0 0.2];

cylEdgeColor = [0.1 0.1 0.1];
cylEdgeAlpha = 0.45;

figure('Name', sprintf('Cube–Cylinder | tol=%.3f', tol), 'Color','w');
ax = gca; hold(ax,'on'); axis(ax,'equal');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
title(ax, sprintf('Cube (master) – Cylinder (slave) | Tol = %.3f mm | A = %.4f mm^2 | Contact tris = %d', ...
                  tol, A_contact, nnz(contactMask)));
view(ax,3); grid(ax,'on');
camlight(ax); lighting(ax,'gouraud');

% cube (master)
patch('Faces', Fc, 'Vertices', Vc, ...
      'FaceColor', cubeColor, ...
      'EdgeColor', 'none', ...
      'FaceAlpha', 0.25, ...
      'Parent', ax);

% cylinder face colouring (mask is SLAVE faces)
nFacesCyl = size(Fs,1);
faceColors = repmat(cylBaseColor, nFacesCyl, 1);
faceColors(contactMask,:) = repmat(contactColor, nnz(contactMask), 1);

patch('Faces', Fs, 'Vertices', Vs_aligned, ...
      'FaceVertexCData', faceColors, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none', ...
      'FaceAlpha', 0.90, ...
      'Parent', ax);

% cylinder mesh edges overlay
%patch('Faces', Fs, 'Vertices', Vs_aligned, ...
 %     'FaceColor', 'none', ...
  %    'EdgeAlpha', cylEdgeAlpha, ...
   %   'LineWidth', 0.5, ...
    %  'Parent', ax);

% ====================== helper ======================
function avgEdge = computeAverageEdgeLength(F, V)
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    v1 = V(E(:,1),:);
    v2 = V(E(:,2),:);
    edgeLens = sqrt(sum((v1 - v2).^2, 2));
    avgEdge = mean(edgeLens);
end
