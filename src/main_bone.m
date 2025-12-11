function main_bone()
% Static femur–tibia contact with adaptive tolerance = baseTol + min vertical gap

    clear; clc; close all;

    scriptDir = fileparts(mfilename('fullpath'));
    if exist(fullfile(scriptDir, 'model'), 'dir')
        projectRoot = scriptDir;
    else
        projectRoot = fileparts(scriptDir);
    end

    modelDir = fullfile(projectRoot, 'model');
    addpath(genpath(fullfile(projectRoot, 'src')));

    fprintf("Project root detected: %s\n", projectRoot);

    % Load knee bones (using MATLAB's stlread, triangulation)
    femurMesh = stlread(fullfile(modelDir, "Femur.STL"));
    tibiaMesh = stlread(fullfile(modelDir, "Tibia.STL"));

    nFrames = 1;

    % 4x4xN transforms, start as identities
    T_femur = repmat(eye(4), 1, 1, nFrames);
    T_tibia = repmat(eye(4), 1, 1, nFrames);

    % Tibia orientation: rotate about global Y (coronal plane)
    y_deg_tibia = 180;  % adjust if CAD tells you otherwise
    Ry = [ cosd(y_deg_tibia)   0   sind(y_deg_tibia);
           0                   1   0;
          -sind(y_deg_tibia)   0   cosd(y_deg_tibia) ];
    T_tibia(1:3,1:3,1) = Ry;

    % Tibia translation: move in Z until condyles are near femur
    % Replace offsetZ_tibia with CAD-measured distance if you have it
    offsetZ_tibia = -50;   % mm example
    T_tibia(1:3,4,1) = [0; 0; offsetZ_tibia];

    % Femur left at identity
    % T_femur(:,:,1) = eye(4);

    % Apply transforms to vertices
    V_femur_frames = applyTransformSeries(femurMesh.Points, T_femur);
    V_tibia_frames = applyTransformSeries(tibiaMesh.Points, T_tibia);

    Vf1 = V_femur_frames(:,:,1);
    Vt1 = V_tibia_frames(:,:,1);

    % Build body structs for contact
    femurBody = buildBodyStruct(femurMesh.ConnectivityList, Vf1);
    tibiaBody = buildBodyStruct(tibiaMesh.ConnectivityList, Vt1);

    % Make sure master has KD-tree (for STS)
    if isfield(femurBody,'triCentroid')
        femurBody.kdtree = KDTreeSearcher(femurBody.triCentroid);
    else
        femurBody.kdtree = KDTreeSearcher(femurBody.V);
    end

    % Contact parameters
    baseTol = 2;   % your chosen tolerance (e.g. cartilage thickness, mm)

    opts.sampleThreshold   = 0.2;
    opts.neighRadiusFactor = 5.0;
    opts.maxNeighbours     = 50;
    opts.roiExpandFactor   = 1.5;
    opts.debugPlot         = false;

    % Adaptive contact: effective tol = baseTol + min vertical gap
    [A_tf, contactMask_tf, tol_eff, e_min] = computeContactArea_adaptiveGapVertical( ...
        tibiaBody, femurBody, baseTol, opts);

    fprintf("Min vertical gap between tibia and femur (Z): %.4f\n", e_min);
    fprintf("Effective tolerance (baseTol + gap): %.4f\n", tol_eff);
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

    % Tibia contact region in neon green
    inContact = contactMask_tf;
    Ft_c = Ft_all(inContact, :);
    if ~isempty(Ft_c)
        patch("Faces", Ft_c, ...
              "Vertices", Vt1, ...
              "FaceColor", [0 1 0], ...
              "EdgeColor", "none", ...
              "FaceAlpha", 0.9);
    end

    axis equal; grid on;
    xlabel("X"); ylabel("Y"); zlabel("Z");
    title(sprintf("Femur–Tibia static pose: A_{TF} = %.3f", A_tf));
    view(60,30); camlight headlight; lighting gouraud;

    % If you want supervisor-style coronal view:
    % view(0,0);

    % Local helper
    function V_frames = applyTransformSeries(V, T)
        nV  = size(V, 1);
        Vh  = [V, ones(nV,1)]';      % 4 x nVerts
        VhN = pagemtimes(T, Vh);     % 4 x nVerts x N
        V_frames = permute(VhN(1:3,:,:), [2 1 3]);  % nVerts x 3 x N
    end
end