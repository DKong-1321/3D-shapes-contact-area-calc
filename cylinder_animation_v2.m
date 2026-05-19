% run_cylinder_visualisation.m
% =====================================================
% Place this file in the root of 3D-shapes-contact-area-calc
% (same level as data/, model/, lib/).
% Run it from that folder in MATLAB.
% =====================================================

%% ---- 1. Add lib/ to path (must happen before load) ----
% The classes (@Digitisation, @JCS, +Knee, etc.) must be on the path
% BEFORE load() is called, otherwise MATLAB cannot reconstruct the
% MCOS objects and throws:
%   "Variable cannot be instantiated as an object"

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'lib')));

assert(exist('Digitisation', 'class') == 8, ...
    'Digitisation class not found — check lib/ is in this project folder');
assert(exist('JCS', 'class') == 8, ...
    'JCS class not found — check lib/ is in this project folder');
fprintf('Classes on path: OK\n');

%% ---- 2. Load matlab_1.mat ----

matFile = fullfile(projectRoot, 'data', 'matlab_1.mat');
assert(exist(matFile, 'file') == 2, 'Cannot find data/matlab_1.mat');

fprintf('Loading matlab_1.mat ...\n');
S = load(matFile);

% Variables inside matlab_1.mat:
%   digitisation  — Digitisation class object
%   kinematics    — Kinematics class object  
%   config        — plain struct
%   module        — Module enum
%   root          — char (original save path, can be ignored)

digitisation = S.digitisation;
kinematics   = S.kinematics;
config       = S.config;
fprintf('Loaded OK.\n');

%% ---- 3. Patch STL paths to this machine ----
% config.stl.femur_right and tibia_right were saved pointing at your
% supervisor's machine. Override to the local model/ folder.

config.stl.femur_right = fullfile(projectRoot, 'model', 'cylinder', 'top_final1.stl');
config.stl.tibia_right = fullfile(projectRoot, 'model', 'cylinder', 'bottom_final1.stl');

assert(exist(config.stl.femur_right, 'file') == 2, ...
    'Cannot find model/cylinder/top_final1.stl');
assert(exist(config.stl.tibia_right, 'file') == 2, ...
    'Cannot find model/cylinder/bottom_final1.stl');

%% ---- 4. Extract transforms for trial 1 ----

TRIAL_INDEX = 1;
tr = kinematics.trajectories(TRIAL_INDEX);

% Transform field names inside each Trajectory object:
%   tr.Transform.gTf  [4x4x201]  global -> femur body frame
%   tr.Transform.gTt  [4x4x201]  global -> tibia body frame
%   tr.Transform.fTt  [4x4x201]  tibia body in femur body frame (~48 mm)

gTf     = tr.Transform.gTf;
gTt     = tr.Transform.gTt;
fTt     = tr.Transform.fTt;
nFrames = size(fTt, 3);

fprintf('Trial %d | %d frames\n', TRIAL_INDEX, nFrames);
fprintf('fTt frame 1 translation: [%.1f  %.1f  %.1f] mm\n', fTt(1:3,4,1)');

% Quick sanity check — fTt should be ~48 mm, not ~1500 mm
if norm(fTt(1:3,4,1)) > 500
    warning('fTt looks like a global transform (%.0f mm) — expected ~48 mm', ...
            norm(fTt(1:3,4,1)));
end

%% ---- 5. Digitisation landmark axes ----
% digitisation.visualise() plots the 3D positions of the medial, lateral
% and distal/proximal digitised points for both cylinders, with
% medial-lateral and proximal-distal axis arrows.

fprintf('\n[1/4] Digitisation landmarks...\n');
digitisation.visualise();

%% ---- 6. Raw tracker frame axes (frame 1) ----
% Checks the femur (Y tracker) and tibia (T tracker) coordinate frames
% look right-handed and point in sensible directions.

fprintf('[2/4] Tracker frame axes at frame 1...\n');
figure(2);
visualise_tracker_transforms(gTf(:,:,1), gTt(:,:,1), config.is_right_knee);
title(sprintf('Tracker coordinate frames — frame 1 (trial %d)', TRIAL_INDEX));

%% ---- 7. Kinematics rotation animation ----
% Animates the femur and tibia rotation vectors across all frames.
% Uses gTf and gTt (global lab-frame transforms).

fprintf('[3/4] Kinematics rotation animation...\n');
transform.femur = gTf;
transform.tibia = gTt;
visualise_kinematics(transform);

%% ---- 8. STL animation ----
% Applies gTf and gTt to the STL vertices and animates.
% NOTE: visualise_stl contains a  while true  loop — it runs
% until you close the figure window. It also saves stlplot.gif.

fprintf('[4/4] STL animation  (close figure window to stop)...\n');
visualise_stl(config, gTt, gTf);