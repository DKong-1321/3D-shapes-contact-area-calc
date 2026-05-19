clear; clc;

%% -------- Paths --------
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;

dataDir  = fullfile(projectRoot,'data');
modelDir = fullfile(projectRoot,'model','cylinder');

matFilePath  = fullfile(dataDir,'cylinders.mat');
femurStlPath = fullfile(modelDir,'top_final1.stl');
tibiaStlPath = fullfile(modelDir,'bottom_final1.stl');

%% -------- Load data --------
S  = load(matFilePath);
tr = S.trajectories(1);

fTt = tr.Transform.fTt;

%% ============================================================
%% 1️⃣ FRAME 1 RELATIVE TRANSFORM CHECK
%% ============================================================

T1 = fTt(:,:,1);

disp('----- fTt(:,:,1) -----');
disp(T1);

% Identity comparison
I4 = eye(4);
identityError = norm(T1 - I4);

fprintf('\nNorm difference from identity: %.6f\n', identityError);

if identityError < 1e-6
    disp('Frame 1 is essentially identity.');
else
    disp('Frame 1 is NOT identity (this may be expected depending on setup).');
end

%% ============================================================
%% 2️⃣ ROTATION ORTHOGONALITY CHECK
%% ============================================================

R = T1(1:3,1:3);
orthError = norm(R'*R - eye(3));

fprintf('\nRotation orthogonality error: %.6e\n', orthError);

if orthError < 1e-10
    disp('Rotation matrix is orthogonal (good rigid transform).');
else
    disp('Rotation matrix is NOT perfectly orthogonal.');
end

detR = det(R);
fprintf('det(R) = %.6f\n', detR);

if abs(detR - 1) < 1e-6
    disp('det(R) ≈ 1 (proper rotation).');
else
    disp('det(R) not equal to 1 (possible scaling or reflection).');
end

%% ============================================================
%% 3️⃣ TRANSLATION MAGNITUDE CHECK
%% ============================================================

t = T1(1:3,4);
tMag = norm(t);

fprintf('\nTranslation vector (frame 1): [%.4f %.4f %.4f]\n', t);
fprintf('Translation magnitude: %.4f\n', tMag);

%% ============================================================
%% 4️⃣ STL SIZE CHECK (UNIT COMPARISON)
%% ============================================================

Sf = stlread(femurStlPath);
Vf = Sf.Points;

bboxMin = min(Vf);
bboxMax = max(Vf);

bboxSize = bboxMax - bboxMin;
modelScale = norm(bboxSize);

fprintf('\nApprox STL bounding box diagonal: %.4f\n', modelScale);

%% ============================================================
%% 5️⃣ UNIT HEURISTIC CHECK
%% ============================================================

if modelScale > 200 && tMag < 1
    disp('⚠ STL likely in mm, transforms likely in meters.');
elseif modelScale < 1 && tMag > 10
    disp('⚠ STL likely in meters, transforms likely in mm.');
else
    disp('Units appear roughly consistent.');
end

disp('----- Diagnostic complete -----');
