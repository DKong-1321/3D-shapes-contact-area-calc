%Paths
dataDir  = fullfile(pwd,'data');
modelDir = fullfile(pwd,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');        % femur (top)
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');   % tibia (bottom)

assert(exist(matPath,'file')==2,  'Missing: %s', matPath);
assert(exist(femurSTL,'file')==2, 'Missing: %s', femurSTL);
assert(exist(tibiaSTL,'file')==2, 'Missing: %s', tibiaSTL);

%Options
trialIndex = 1;
useFlexionOffset  = true;
useAnteriorOffset = true;

flipFlexionSign  = false;
flipAnteriorSign = true;

framePause = 0.02;

%MAT
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

%First row of kinematics
K = tr.Kinematics;
if isfield(K,'tibiofemoral') && ~isempty(K.tibiofemoral)
    TF = K.tibiofemoral;
elseif isfield(K,'tibiafemoral') && ~isempty(K.tibiafemoral)
    TF = K.tibiafemoral;
else
    error('No tibiofemoral kinematics found in trajectories(%d).Kinematics', trialIndex);
end

if ~istable(TF)
    error('Expected tibiofemoral kinematics to be a table.');
end

row1 = TF(1,:);
vars = TF.Properties.VariableNames;

flexion  = row1{1, strcmpi(vars,'flexion')};    % deg
anterior = row1{1, strcmpi(vars,'anterior')};   % mm

if flipFlexionSign,  flexion  = -flexion;  end
if flipAnteriorSign, anterior = -anterior; end

fprintf('Row1 flexion  = %.4f deg\n', flexion);
fprintf('Row1 anterior = %.4f mm\n',  anterior);

%Constant offset about origin
Tcorr = eye(4);

if useFlexionOffset
    th = deg2rad(flexion);
    R = [ 1  0        0;
          0  cos(th) -sin(th);
          0  sin(th)  cos(th) ];
    Tcorr(1:3,1:3) = R;
end

if useAnteriorOffset
    Tcorr(1:3,4) = [0; anterior; 0];
end

TcorrInv = inv(Tcorr);

%Motion matrices (GLOBAL)
Tstruct = tr.Transform;

if ~isfield(Tstruct,'gTfi') || isempty(Tstruct.gTfi)
    error('No gTfi found in trajectories(%d).Transform', trialIndex);
end
if ~isfield(Tstruct,'gTti') || isempty(Tstruct.gTti)
    error('No gTti found in trajectories(%d).Transform', trialIndex);
end

gTfi = Tstruct.gTfi;   % femur in global
gTti = Tstruct.gTti;   % tibia in global

nFrames = size(gTfi,3);

%APPLY OFFSET TO THE MATRICES
gTfi_corr = zeros(size(gTfi));
for k = 1:nFrames
    gTfi_corr(:,:,k) = gTfi(:,:,k) * TcorrInv;   % post-multiply correction
end

%LOAD STLs
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm = Sf.ConnectivityList;
Vt0 = St.Points;  Ft = St.ConnectivityList;

%PLOT SETUP
figure('Name','Global (camera) frame, both bodies moving (offset applied)'); hold on;
axis equal; grid on; view(3);
xlabel('X (medial)'); ylabel('Y (anterior)'); zlabel('Z (superior)');

hTibia = patch('Faces',Ft,'Vertices',Vt0,'FaceAlpha',0.20,'EdgeColor','none');
hFemur = patch('Faces',Fm,'Vertices',Vm0,'FaceAlpha',0.50,'EdgeColor','none');

%ANIMATE
for k = 1:nFrames
    % tibia in global
    Tt = gTti(:,:,k);
    Rt = Tt(1:3,1:3);
    tt = Tt(1:3,4)';

    Vt_k = (Rt * Vt0')' + tt;

    % femur in global (corrected)
    Tf = gTfi_corr(:,:,k);
    Rf = Tf(1:3,1:3);
    tf = Tf(1:3,4)';

    Vm_k = (Rf * Vm0')' + tf;

    set(hTibia,'Vertices',Vt_k);
    set(hFemur,'Vertices',Vm_k);

    title(sprintf('Trial %d | Frame %d / %d', trialIndex, k, nFrames));
    drawnow;
    pause(framePause);
end
