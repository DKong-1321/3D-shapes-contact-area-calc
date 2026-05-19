%% GLOBAL CONTACT + CONTACT AREA HISTORY
clear; clc; close all;

%% ================= PATHS =================
dataDir  = fullfile(pwd,'data');
modelDir = fullfile(pwd,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');        % femur
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');   % tibia

assert(exist(matPath,'file')==2,  'Missing: %s', matPath);
assert(exist(femurSTL,'file')==2, 'Missing: %s', femurSTL);
assert(exist(tibiaSTL,'file')==2, 'Missing: %s', tibiaSTL);

%% ================= OPTIONS =================
trialIndex = 1;
tol        = 0.2;        % contact tolerance [mm]
framePause = 0.02;

useFlexionOffset  = true;
useAnteriorOffset = true;

flipFlexionSign  = false;
flipAnteriorSign = true;

%% ================= LOAD MAT =================
S = load(matPath);

trajectories = [];
if isfield(S,'trajectories'), trajectories = S.trajectories; end
if isempty(trajectories) && isfield(S,'Trajectories'), trajectories = S.Trajectories; end
if isempty(trajectories)
    fns = fieldnames(S);
    trajectories = S.(fns{1});
end
if istable(trajectories), trajectories = table2array(trajectories); end
if iscell(trajectories),  trajectories = [trajectories{:}]; end

tr = trajectories(trialIndex);

%% ================= KINEMATICS OFFSET =================
K = tr.Kinematics;
if isfield(K,'tibiofemoral') && ~isempty(K.tibiofemoral)
    TF = K.tibiofemoral;
elseif isfield(K,'tibiafemoral') && ~isempty(K.tibiafemoral)
    TF = K.tibiafemoral;
else
    error('No tibiofemoral kinematics found');
end

row1 = TF(1,:);
vars = TF.Properties.VariableNames;

flexion  = row1{1, strcmpi(vars,'flexion')};    % deg
anterior = row1{1, strcmpi(vars,'anterior')};   % mm

if flipFlexionSign,  flexion  = -flexion;  end
if flipAnteriorSign, anterior = -anterior; end

fprintf('Initial flexion  = %.3f deg\n', flexion);
fprintf('Initial anterior = %.3f mm\n', anterior);

%% ================= OFFSET TRANSFORM =================
Tcorr = eye(4);

if useFlexionOffset
    th = deg2rad(flexion);
    Tcorr(1:3,1:3) = [ 1  0        0;
                       0  cos(th) -sin(th);
                       0  sin(th)  cos(th) ];
end

if useAnteriorOffset
    Tcorr(1:3,4) = [0; anterior; 0];
end

TcorrInv = inv(Tcorr);

%% ================= GLOBAL MOTION MATRICES =================
Tstruct = tr.Transform;

assert(isfield(Tstruct,'gTfi') && ~isempty(Tstruct.gTfi),'Missing gTfi');
assert(isfield(Tstruct,'gTti') && ~isempty(Tstruct.gTti),'Missing gTti');

gTfi = Tstruct.gTfi;   % femur in global
gTti = Tstruct.gTti;   % tibia in global

nFrames = size(gTfi,3);

% Apply offset to femur matrices only
gTfi_corr = zeros(size(gTfi));
for k = 1:nFrames
    gTfi_corr(:,:,k) = gTfi(:,:,k) * TcorrInv;
end

%% ================= LOAD STLs =================
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm = Sf.ConnectivityList;
Vt0 = St.Points;  Ft = St.ConnectivityList;

%% ================= MASTER BODY (TIBIA) =================
master = buildBodyStruct(Ft, Vt0);

%% ================= VISUALISATION SETUP =================
figure('Name','Global contact animation'); hold on;
axis equal; grid on; view(3);
xlabel('X (medial)'); ylabel('Y (anterior)'); zlabel('Z (superior)');

hTibia = patch('Faces',Ft,'Vertices',Vt0,...
    'FaceAlpha',0.20,'EdgeColor','none');

hFemur = patch('Faces',Fm,'Vertices',Vm0,...
    'FaceColor','flat','EdgeColor','none');

%% ================= STORAGE =================
contactArea = zeros(nFrames,1);

%% ================= ANIMATION LOOP =================
for k = 1:nFrames

    % ---- Tibia (global) ----
    Tt = gTti(:,:,k);
    Vt = (Tt(1:3,1:3) * Vt0')' + Tt(1:3,4)';
    set(hTibia,'Vertices',Vt);

    % ---- Femur (global, corrected) ----
    Tf = gTfi_corr(:,:,k);
    Vm = (Tf(1:3,1:3) * Vm0')' + Tf(1:3,4)';
    set(hFemur,'Vertices',Vm);

    % ---- Slave body (femur) ----
    slave = buildBodyStruct(Fm, Vm);

    % ---- Contact ----
    res = computeContactArea_STS_hybrid(master, slave, tol);
    contactArea(k) = res.contactArea;

    % ---- Colour contact triangles ----
    C = repmat([0.6 0.6 0.6], size(Fm,1), 1);
    if any(res.contactMask)
        C(res.contactMask,:) = repmat([1 0 0], nnz(res.contactMask), 1);
    end
    set(hFemur,'FaceVertexCData',C);

    title(sprintf('Frame %d / %d | Contact = %.2f mm^2',...
        k, nFrames, contactArea(k)));

    drawnow;
    pause(framePause);
end

%% ================= CONTACT AREA HISTORY =================
figure;
plot(contactArea,'LineWidth',2);
grid on;
xlabel('Frame');
ylabel('Contact area (mm^2)');
title(sprintf('Contact area vs time (trial %d)', trialIndex));

%% ================= HELPER =================
function body = buildBodyStruct(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);

    body.F = F;
    body.V = V;
    body.centroids = (v1 + v2 + v3) / 3;
    body.areas = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
end
