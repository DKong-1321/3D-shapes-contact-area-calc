% cylinder_animation.m
% Standalone STL animation for the cylinder experiment.
% Run from the root of 3D-shapes-contact-area-calc.
% Requires: data/matlab 1.mat  (must contain kinematics.trajectories)
%           model/cylinder/top.stl    (femur — top cylinder)
%           model/cylinder/bottom.stl (tibia — bottom cylinder)
%
% The STL files have their origin at a bounding-box corner, not at the
% body frame origin used in the optical tracking pipeline.
% Body frame origin = midpoint of medial and lateral probe touches = centre
% of the mating face of each cylinder.
% This script recentres each mesh onto its mating face centre before
% applying the optical tracking transforms, so the two cylinders sit
% together correctly at frame 1.
%
% top.stl    mating face (bottom of top cylinder) is at Z_min ~ 0.42 mm
%            XY centre of that face ~ (20.28, 38.51)
%            Correction applied: subtract (20.28, 38.51, 0.42) from all vertices
%
% bottom.stl mating face (top of bottom cylinder) is at Z_max ~ 83.50 mm
%            XY centre of that face ~ (20.28, 37.22)
%            Correction applied: subtract (20.28, 37.22, 83.50) from all vertices
%
% After correction both mesh origins sit at their respective mating face centres,
% matching the body frame origin convention from digitisation.

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'lib')));
addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'external')));

% load transforms
S          = load(fullfile(projectRoot, 'data', 'matlab 1.mat'));
kinematics = S.kinematics;
TRIAL      = 1;
tr         = kinematics.trajectories(TRIAL);
gTfi       = tr.Transform.gTfi;   % [4x4xN] femur body in global
gTti       = tr.Transform.gTti;   % [4x4xN] tibia body in global
fTt        = tr.Transform.fTt;    % [4x4xN] tibia in femur frame
nFrames    = size(fTt, 3);

fprintf('Trial %d | %d frames\n', TRIAL, nFrames);
fprintf('fTt frame 1 translation norm: %.1f mm\n', norm(fTt(1:3,4,1)));

% compute tTf = fTt^-1 (femur in tibia frame) using rigid inversion
tTf = zeros(4, 4, nFrames);
for k = 1:nFrames
    R = fTt(1:3, 1:3, k);
    t = fTt(1:3, 4,   k);
    tTf(:,:,k) = [R' -R'*t; 0 0 0 1];
end

% load STLs
femur_mesh = stlread(fullfile(projectRoot, 'model', 'cylinder', 'top.stl'));
tibia_mesh = stlread(fullfile(projectRoot, 'model', 'cylinder', 'bottom.stl'));

V_femur = femur_mesh.Points;
F_femur = femur_mesh.ConnectivityList;
V_tibia = tibia_mesh.Points;
F_tibia = tibia_mesh.ConnectivityList;

% recentre each mesh onto its mating face centre
% top (femur): mating face is the bottom face, Z_min side
femur_origin = [mean(V_femur(:,1)),  mean(V_femur(:,2)),  min(V_femur(:,3))];
V_femur = V_femur - femur_origin;

% bottom (tibia): mating face is the top face, Z_max side
tibia_origin = [mean(V_tibia(:,1)),  mean(V_tibia(:,2)),  max(V_tibia(:,3))];
V_tibia = V_tibia - tibia_origin;

% apply rigid transform helper
apply = @(T, V) (T(1:3,1:3) * V')' + T(1:3,4)';

% figure setup
fig = figure('Color', 'w', 'Name', 'Cylinder animation');
axis equal; grid on; view(30, 20);
xlabel('X'); ylabel('Y'); zlabel('Z');
camlight; lighting gouraud;

% draw both meshes at frame 1 to set axis limits
V_femur_0 = apply(tTf(:,:,1), V_femur);
hFemur = patch('Faces', F_femur, 'Vertices', V_femur_0, ...
    'FaceColor', [0.85 0.72 0.60], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hTibia = patch('Faces', F_tibia, 'Vertices', V_tibia, ...
    'FaceColor', [0.60 0.75 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
legend([hFemur hTibia], {'Femur (top)', 'Tibia (bottom)'});
title(sprintf('Frame 1 / %d', nFrames));

% check separation at frame 1
sep = norm(tTf(1:3,4,1));
fprintf('Frame 1 femur origin in tibia frame: [%.2f  %.2f  %.2f] mm  (norm = %.1f mm)\n', ...
    tTf(1,4,1), tTf(2,4,1), tTf(3,4,1), sep);
if sep > 100
    warning('Large separation at frame 1 (%.1f mm) — STL origin correction may need adjusting.', sep);
end

% animate
framePause = 0.03;
gifFile    = 'cylinder_animation.gif';
saveGif    = true;

for k = 1:nFrames
    V_frame = apply(tTf(:,:,k), V_femur);
    set(hFemur, 'Vertices', V_frame);
    title(sprintf('Frame %d / %d', k, nFrames));
    drawnow;

    if saveGif
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, 'DelayTime', framePause);
        else
            imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', framePause);
        end
    end

    pause(framePause);
end

fprintf('Done. GIF saved to %s\n', gifFile);