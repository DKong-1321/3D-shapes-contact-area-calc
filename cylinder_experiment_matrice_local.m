%% LOCAL CONTACT ANIMATION (fTt, tibia fixed)
clear; clc; close all;

%Paths
dataDir  = fullfile(pwd,'data');
modelDir = fullfile(pwd,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');

assert(exist(matPath,'file')==2);
assert(exist(femurSTL,'file')==2);
assert(exist(tibiaSTL,'file')==2);

%Options
trialIndex = 1;
tol = 0.2;
framePause = 0.02;

%% Load MAT
S = load(matPath);
trajectories = S.trajectories;
tr = trajectories(trialIndex);

%% Kinematics (first row offset)
TF = tr.Kinematics.tibiofemoral;
flexion  = TF.flexion(1);
anterior = TF.anterior(1);

Tcorr = eye(4);
th = deg2rad(flexion);
Tcorr(1:3,1:3) = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
Tcorr(1:3,4)   = [0; anterior; 0];
TcorrInv = inv(Tcorr);

%% Transforms
fTt = tr.Transform.fTt;   % tibia in femur
nFrames = size(fTt,3);

%% Load STLs
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm = Sf.ConnectivityList;
Vt  = St.Points;  Ft = St.ConnectivityList;

%% Build MASTER once (tibia)
master = buildBodyStruct(Ft, Vt);

%% Plot
figure; hold on; axis equal; grid on; view(3);
patch('Faces',Ft,'Vertices',Vt,'FaceAlpha',0.2,'EdgeColor','none');

hFemur = patch('Faces',Fm,'Vertices',Vm0,...
               'FaceColor','flat','EdgeColor','none');

%% Animate
for k = 1:nFrames
    % Correct transform
    fTt_corr = TcorrInv * fTt(:,:,k);

    % femur in tibia CS
    tTf = inv(fTt_corr);

    R = tTf(1:3,1:3);
    t = tTf(1:3,4)';

    Vm = (R * Vm0')' + t;
    set(hFemur,'Vertices',Vm);

    % SLAVE body (rebuild per frame)
    slave = buildBodyStruct(Fm, Vm);

    % Contact
    res = computeContactArea_STS_hybrid(master, slave, tol);

    % Colour faces
    C = repmat([0.7 0.7 0.7], size(Fm,1), 1);
    C(res.contactMask,:) = repmat([1 0 0], nnz(res.contactMask), 1);
    set(hFemur,'FaceVertexCData',C);

    title(sprintf('Frame %d | A = %.2f mm^2', k, res.contactArea));
    drawnow;
    pause(framePause);
end

%% ---- helper ----
function body = buildBodyStruct(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    body.F = F;
    body.V = V;
    body.centroids = (v1+v2+v3)/3;
    body.areas = 0.5*vecnorm(cross(v2-v1,v3-v1,2),2,2);
end
