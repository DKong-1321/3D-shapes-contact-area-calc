% main_cube_cylinder.m
% Cube–cylinder contact for two positions over a range of tolerances

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);  % go up one level
end

addpath(genpath(fullfile(projectRoot, 'src')));  % helpers, transforms, etc.

modelDir = fullfile(projectRoot, 'model');

cubeFile = fullfile(modelDir, '50mm cube 4 iterations.stl');
cylFile  = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');

tolList = [0.05 0.10 0.20 0.30 0.40 0.50];  % tolerances to test (mm)

opts.sampleThreshold   = 0.5;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

% Load meshes
[Fc, Vc] = loadStlMesh(cubeFile, 'cube');
[Fs, Vs] = loadStlMesh(cylFile,  'cylinder');

% Place cylinder on cube once (reference configuration)
[Vs_aligned, cylShift] = placeCylinderOnCube(Vs, Vc);

cubeRef_V     = Vc;
cylinderRef_V = Vs_aligned;

cubeRef     = buildBodyStruct(Fc, cubeRef_V);
cylinderRef = buildBodyStruct(Fs, cylinderRef_V);

fprintf('Cube:     %d faces, %d verts\n', size(Fc,1), size(cubeRef_V,1));
fprintf('Cylinder: %d faces, %d verts\n', size(Fs,1), size(cylinderRef_V,1));

% Mesh resolutions (average edge length)
cubeEdge = computeAverageEdgeLength(Fc, cubeRef_V);
cylEdge  = computeAverageEdgeLength(Fs, cylinderRef_V);
fprintf('Average edge length (cube):     %.4f mm\n', cubeEdge);
fprintf('Average edge length (cylinder): %.4f mm\n', cylEdge);

% Two positions for the cylinder: base and shifted in +X
Npos = 2;
T_cube = repmat(eye(4), 1, 1, Npos);  % cube fixed

T_cyl = repmat(eye(4), 1, 1, Npos);
T_cyl(:,:,1) = eye(4);                % Position 1: base pose

T2 = eye(4);
T2(1,4) = 20.0;                       % Position 2: shift +20 mm in X
T_cyl(:,:,2) = T2;

posLabels = { ...
    'Position 1: base (central contact)', ...
    'Position 2: shifted +20 mm in X'};

% Colours
cubeColor    = [0.8 0.8 0.8];
cylBaseColor = [0.25 0.25 0.25];    % dark grey
contactColor = [0.0 1.0 0.0];       % neon green

for it = 1:numel(tolList)
    tol = tolList(it);
    fprintf('\n===== Tolerance = %.3f mm =====\n', tol);

    for p = 1:Npos
        % Transform vertices for this position
        Vc_k = transformVertices(cubeRef_V,     T_cube(:,:,p));
        Vs_k = transformVertices(cylinderRef_V, T_cyl(:,:,p));

        % Build body structs
        cube_k     = buildBodyStruct(Fc, Vc_k);
        cylinder_k = buildBodyStruct(Fs, Vs_k);

        master = cube_k;      % cube = master
        slave  = cylinder_k;  % cylinder = slave

        master.kdtree = KDTreeSearcher(master.triCentroid);

        % Contact calculation
        [A_contact, contactMask] = computeContactArea_STS(slave, master, tol, opts);

        fprintf('%s: contact area = %.6f mm^2, contact tris = %d\n', ...
            posLabels{p}, A_contact, nnz(contactMask));

        % Per-tolerance, per-position figure
        figName = sprintf('tol = %.3f, %s', tol, posLabels{p});
        figure('Name', figName, 'Color','w');
        ax = gca; hold(ax,'on');
        axis(ax,'equal');
        xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
        title(ax, sprintf('%s, tol = %.3f, A = %.4f mm^2', ...
                          posLabels{p}, tol, A_contact));
        view(ax,3); camlight(ax); lighting(ax,'gouraud'); grid(ax,'on');

        patch('Faces', Fc, 'Vertices', Vc_k, ...
              'FaceColor', cubeColor, ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.3, ...
              'Parent', ax);

        nFacesSlave = size(slave.F,1);
        faceColors  = repmat(cylBaseColor, nFacesSlave, 1);
        faceColors(contactMask,:) = repmat(contactColor, nnz(contactMask), 1);

        patch('Faces', Fs, 'Vertices', Vs_k, ...
              'FaceVertexCData', faceColors, ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.95, ...
              'Parent', ax);
    end
end

% Helper: apply 4x4 transform to vertices
function V_out = transformVertices(V_in, T)
    Vh = [V_in, ones(size(V_in,1),1)];
    Vw = (T * Vh.').';
    V_out = Vw(:,1:3);
end

% Helper: average edge length of a triangular mesh
function avgEdge = computeAverageEdgeLength(F, V)
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    E = unique(E,'rows');
    v1 = V(E(:,1),:);
    v2 = V(E(:,2),:);
    edgeVecs = v1 - v2;
    edgeLens = sqrt(sum(edgeVecs.^2, 2));
    avgEdge = mean(edgeLens);
end