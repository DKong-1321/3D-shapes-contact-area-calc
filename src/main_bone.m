<<<<<<< Updated upstream
function main_bone()
% Use optical transform data as ground-truth poses for femur, tibia,
% and patella, optionally re-orient the whole knee, then run contact.
=======
% main_bone.m
% Static femur–tibia contact with adaptive tolerance based on min surface gap.
>>>>>>> Stashed changes

    clear; clc;
    scriptDir = fileparts(mfilename('fullpath'));

<<<<<<< Updated upstream
    if exist(fullfile(scriptDir, 'model'), 'dir')
        projectRoot = scriptDir;
    else
        projectRoot = fileparts(scriptDir);
    end

    modelDir = fullfile(projectRoot, 'model');
    dataDir  = fullfile(projectRoot, 'data');

    addpath(genpath(fullfile(projectRoot, 'src')));

    fprintf('Loading STL meshes...\n');
    [F_fem, V_fem] = loadAnyStl(fullfile(modelDir, 'Femur.STL'));
    [F_tib, V_tib] = loadAnyStl(fullfile(modelDir, 'Tibia.STL'));
    [F_pat, V_pat] = loadAnyStl(fullfile(modelDir, 'Patella.STL'));

    % Load optical transforms
    tfFile = fullfile(dataDir, 'transforms_optical-data.mat');
    fprintf('Loading transform data from %s\n', tfFile);
    S = load(tfFile);

    % Assuming structure: transforms.Intact.in_global.(femur/tibia/patella)
    Tf = S.transforms.Intact.in_global.femur;
    Tt = S.transforms.Intact.in_global.tibia;
    Tp = S.transforms.Intact.in_global.patella;

    nFrames = min([size(Tf,3), size(Tt,3), size(Tp,3)]);
    fprintf('Found %d frames of transform data.\n', nFrames);

    % Choose frame index
    frameIdx = round(nFrames/2);
    fprintf('Using frame %d as ground truth.\n', frameIdx);

    T_fem_GT = Tf(:,:,frameIdx);
    T_tib_GT = Tt(:,:,frameIdx);
    T_pat_GT = Tp(:,:,frameIdx);

    % Apply transforms to meshes
    V_fem_world = applySingleTransform(V_fem, T_fem_GT);
    V_tib_world = applySingleTransform(V_tib, T_tib_GT);
    V_pat_world = applySingleTransform(V_pat, T_pat_GT);

    % Re-orientation using tibia
    useReorientation = true;   % set false to stay in optical global

    if useReorientation
        [V_fem_plot, V_tib_plot, V_pat_plot] = ...
            reorientKneeWithTibia(V_fem_world, V_tib_world, V_pat_world);
        coordLabel = sprintf('optical frame %d (re-oriented)', frameIdx);
    else
        V_fem_plot = V_fem_world;
        V_tib_plot = V_tib_world;
        V_pat_plot = V_pat_world;
        coordLabel = sprintf('optical frame %d (raw global)', frameIdx);
    end

    % Visualise posed knee
    figure; hold on; axis equal; grid on;
    patch('Faces', F_tib, 'Vertices', V_tib_plot, ...
          'FaceColor',[0.8 0.8 0.9],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', F_fem, 'Vertices', V_fem_plot, ...
          'FaceColor',[0.9 0.7 0.7],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', F_pat, 'Vertices', V_pat_plot, ...
          'FaceColor',[0.7 0.9 0.7],'EdgeColor','none','FaceAlpha',0.7);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Knee at ' coordLabel]);
    view(3); camlight; lighting gouraud; rotate3d on;

    % Build body structs and compute contact
    femurBody   = buildBodyStruct(F_fem, V_fem_plot,   'femur');
    tibiaBody   = buildBodyStruct(F_tib, V_tib_plot,   'tibia');
    patellaBody = buildBodyStruct(F_pat, V_pat_plot,   'patella');

    % Make sure masters have KD-trees
    tibiaBody   = addKDTreeToBody(tibiaBody);
    patellaBody = addKDTreeToBody(patellaBody);

    tol  = 0.25;   % in STL units
    opts.sampleThreshold   = 0.2;
    opts.neighRadiusFactor = 5.0;
    opts.maxNeighbours     = 50;

    fprintf('Computing contact areas (femur–tibia and femur–patella)...\n');
    [A_FT, mask_FT] = computeContactArea_STS(femurBody, tibiaBody,   tol, opts); %#ok<NASGU>
    [A_FP, mask_FP] = computeContactArea_STS(femurBody, patellaBody, tol, opts); %#ok<NASGU>

    fprintf('Femur–Tibia  contact area: %.3f\n', A_FT);
    fprintf('Femur–Patella contact area: %.3f\n', A_FP);
end


function V_out = applySingleTransform(V_in, T)
    N  = size(V_in, 1);
    Vh = [V_in, ones(N,1)];      % N x 4
    Vt = (T * Vh.').';           % N x 4
    V_out = Vt(:,1:3);
end

function [V_fem_out, V_tib_out, V_pat_out] = ...
         reorientKneeWithTibia(V_fem, V_tib, V_pat)

    allVerts = [V_fem; V_tib; V_pat];
    centAll  = mean(allVerts, 1);

    V_fem = V_fem - centAll;
    V_tib = V_tib - centAll;
    V_pat = V_pat - centAll;

    C        = cov(V_tib);
    [E, D]   = eig(C);
    [~, idx] = sort(diag(D), 'descend');
    E        = E(:, idx);

    shaftAxis = E(:,1);
    other1    = E(:,2);
    other2    = E(:,3);

    proj1 = V_tib * other1;
    proj2 = V_tib * other2;
    if (max(proj1) - min(proj1)) >= (max(proj2) - min(proj2))
        yAxis = other1;
        xAxis = cross(yAxis, shaftAxis);
    else
        yAxis = other2;
        xAxis = cross(yAxis, shaftAxis);
    end

    xAxis = xAxis / norm(xAxis);
    yAxis = yAxis / norm(yAxis);
    zAxis = shaftAxis / norm(shaftAxis);

    R = [xAxis, yAxis, zAxis];
    if det(R) < 0
        xAxis = -xAxis;
        R     = [xAxis, yAxis, zAxis];
    end

    V_fem = V_fem * R;
    V_tib = V_tib * R;
    V_pat = V_pat * R;

    tibZ     = V_tib(:,3);
    zPlateau = max(tibZ);

    V_fem(:,3) = V_fem(:,3) - zPlateau;
    V_tib(:,3) = V_tib(:,3) - zPlateau;
    V_pat(:,3) = V_pat(:,3) - zPlateau;

    V_fem_out = V_fem;
    V_tib_out = V_tib;
    V_pat_out = V_pat;
end

function [F, V] = loadAnyStl(filename)
    if exist('loadStlMesh','file') == 2
        [F, V] = loadStlMesh(filename);
    else
        data = stlread(filename);
        if isa(data,'triangulation')
            F = data.ConnectivityList;
            V = data.Points;
        else
            F = data.faces;
            V = data.vertices;
        end
    end
end

function body = addKDTreeToBody(body)
    if isfield(body, 'kdtree')
        return;
    end

    F = body.F;
    V = body.V;
    C = ( V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:) ) / 3;
    body.kdtree = KDTreeSearcher(C);
=======
scriptDir = fileparts(mfilename('fullpath'));
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end

modelDir = fullfile(projectRoot, 'model');
addpath(genpath(fullfile(projectRoot, 'src')));

fprintf("Project root detected: %s\n", projectRoot);

% Load knee bones (use the exact filenames you have)
femurMesh = stlread(fullfile(modelDir, "Femur.STL"));
tibiaMesh = stlread(fullfile(modelDir, "Tibia.STL"));

nFrames = 1;

% 4x4xN transforms, start as identities
T_femur = repmat(eye(4), 1, 1, nFrames);
T_tibia = repmat(eye(4), 1, 1, nFrames);

% --- Tibia orientation: rotate in coronal plane about Y ---
y_deg_tibia = 180;  % adjust if CAD tells you otherwise
Ry = [ cosd(y_deg_tibia)   0   sind(y_deg_tibia);
       0                   1   0;
      -sind(y_deg_tibia)   0   cosd(y_deg_tibia) ];
T_tibia(1:3,1:3,1) = Ry;

% --- Tibia translation: move in Z until condyles are near femur ---
% Replace this with the CAD-measured value: origin->surface + cartilage.
offsetZ_tibia = -50;   % mm example; tune this from CAD
T_tibia(1:3,4,1) = [0; 0; offsetZ_tibia];

% Femur left at identity for now
% T_femur(:,:,1) = eye(4);

% Apply transforms to vertices
V_femur_frames = applyTransformSeries(femurMesh.Points, T_femur);
V_tibia_frames = applyTransformSeries(tibiaMesh.Points, T_tibia);

Vf1 = V_femur_frames(:,:,1);
Vt1 = V_tibia_frames(:,:,1);

% Build body structs for contact
femurBody = buildBodyStruct(femurMesh.ConnectivityList, Vf1);
tibiaBody = buildBodyStruct(tibiaMesh.ConnectivityList, Vt1);

% Contact parameters
baseTol = 4;

opts.sampleThreshold   = 0.2;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 1000;
opts.debugPlot         = false;

% Adaptive contact: min gap + baseTol
%[A_tf, contactMask_tf, tol_eff, d_min] = computeContactArea_adaptiveGap( ...
 %   tibiaBody, femurBody, baseTol, opts);

%fprintf("Min 3D gap between tibia and femur: %.4f\n", d_min);
%fprintf("Effective tolerance (gap + baseTol): %.4f\n", tol_eff);
%fprintf("Tibia–Femur contact area (static pose): %.3f\n", A_tf);

% Adaptive contact: min vertical gap + baseTol
[A_tf, contactMask_tf, tol_eff, e_min] = computeContactArea_adaptiveGapVertical( ...
    tibiaBody, femurBody, baseTol, opts);

fprintf("Min vertical gap between tibia and femur (Z): %.4f\n", e_min);
fprintf("Effective tolerance (gap + baseTol): %.4f\n", tol_eff);
fprintf("Tibia–Femur contact area (static pose): %.3f\n", A_tf);


% Visualise geometry and contact
figure("Name","Static femur–tibia alignment","Color","w"); hold on;

% Femur
patch("Faces", femurMesh.ConnectivityList, ...
      "Vertices", Vf1, ...
      "FaceColor", [0.8 0.8 0.8], ...
      "EdgeColor", "none", ...
      "FaceAlpha", 0.8);

% Tibia (base color)
Ft_all = tibiaMesh.ConnectivityList;
patch("Faces", Ft_all, ...
      "Vertices", Vt1, ...
      "FaceColor", [0.2 0.4 0.9], ...
      "EdgeColor", "none", ...
      "FaceAlpha", 0.4);

% Tibia contact region in neon-ish green
inContact = contactMask_tf;
Ft_c = Ft_all(inContact, :);
if ~isempty(Ft_c)
    patch("Faces", Ft_c, ...
          "Vertices", Vt1, ...
          "FaceColor", [0 1 0], ...  % bright green
          "EdgeColor", "none", ...
          "FaceAlpha", 0.9);
end

axis equal; grid on;
xlabel("X"); ylabel("Y"); zlabel("Z");
title(sprintf("Femur–Tibia static pose: A_{TF} = %.3f", A_tf));
view(60,30); camlight headlight; lighting gouraud;

% If you want the supervisor's coronal-style view:
% view(0,0);

% Local helper
function V_frames = applyTransformSeries(V, T)
    nV  = size(V, 1);
    Vh  = [V, ones(nV,1)]';      % 4 x nVerts
    VhN = pagemtimes(T, Vh);     % 4 x nVerts x N
    V_frames = permute(VhN(1:3,:,:), [2 1 3]);  % nVerts x 3 x N
>>>>>>> Stashed changes
end