clear; clc; close all;

%Paths
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;
addpath(genpath(fullfile(projectRoot,'src')));

dataDir  = fullfile(projectRoot,'data');
modelDir = fullfile(projectRoot,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'top_final1.STL');
tibiaSTL = fullfile(modelDir,'bottom_final1.STL');

assert(exist(matPath,'file')==2);
assert(exist(femurSTL,'file')==2);
assert(exist(tibiaSTL,'file')==2);


trialIndex  = 1;
tol         = 0.20;  % mm
framePause  = 0.03;
doDebugPlot = true;

S  = load(matPath);
tr = S.trajectories(trialIndex);

fTt     = tr.Transform.fTt;
nFrames = size(fTt,3);

Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm = Sf.ConnectivityList;
Vt  = St.Points;  Ft = St.ConnectivityList;

master = buildBodyStruct(Ft, Vt);
master.faceNormals = faceNormalsFromFV(master.F, master.V);
master.kdtreeFaces = KDTreeSearcher(master.centroids);
master.bbox        = bboxFromV(master.V);


flex0 = getFirstFrameFlexion(tr);
Rcorr = rotxd(-flex0);
fprintf('Initial flexion offset removed: %.4f deg\n', flex0);


tTf1 = inv(fTt(:,:,1));
tTf1(1:3,1:3) = Rcorr * tTf1(1:3,1:3);

Tstl = stlFrameFromPCA(Vm0);
Treg = tTf1 / Tstl;


if doDebugPlot
    figure('Color','w'); hold on; axis equal; grid on; view(3);
    title('Frame 1 alignment check');
    patch('Faces',Ft,'Vertices',Vt,'FaceAlpha',0.25,'EdgeColor','none');
    Vm_dbg = applyTformToVertices(Treg * tTf1, Vm0);
    patch('Faces',Fm,'Vertices',Vm_dbg,'FaceAlpha',0.5,'EdgeColor','none');
    camlight; lighting gouraud;
end


A = zeros(nFrames,1);

figure('Color','w'); hold on; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');

patch('Faces',Ft,'Vertices',Vt, ...
    'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeColor','none');

hFemur = patch('Faces',Fm,'Vertices',Vm0, ...
    'FaceColor','flat','EdgeColor','none');

colormap([0.7 0.7 0.7; 1 0 1]);
clim([1 2]);

for k = 1:nFrames
    tTf = inv(fTt(:,:,k));
    tTf(1:3,1:3) = Rcorr * tTf(1:3,1:3);
    T  = Treg * tTf;

    Vm = applyTformToVertices(T, Vm0);
    slave = buildBodyStruct(Fm, Vm);

    res = computeContactArea_STS_hybrid(master, slave, tol);
    A(k) = res.contactArea;

    set(hFemur,'Vertices',Vm);

    C = ones(size(Fm,1),1);
    if isfield(res,'contactMask')
        C(res.contactMask) = 2;
    end
    set(hFemur,'FaceVertexCData',C);

    title(sprintf('Frame %d / %d | Contact area = %.3f mm^2', ...
        k, nFrames, A(k)));

    drawnow;
    pause(framePause);
end


figure('Color','w');
plot(A,'-o','LineWidth',1.5);
xlabel('Frame');
ylabel('Contact area (mm^2)');
title('Cylinder contact area vs frame (hybrid)');
grid on;

fprintf('\nMean contact area = %.4f mm^2\n', mean(A));



function flex0 = getFirstFrameFlexion(tr)
TF = tr.Kinematics.tibiofemoral;
if istable(TF)
    flex0 = TF{1,'flexion'};
else
    flex0 = TF(1,1);
end
end

function body = buildBodyStruct(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
body.F = F;
body.V = V;
body.centroids = (v1+v2+v3)/3;
body.areas = 0.5 * vecnorm(cross(v2-v1,v3-v1,2),2,2);
end

function Vout = applyTformToVertices(T,V)
Vh = [V ones(size(V,1),1)];
Vt = (T * Vh')';
Vout = Vt(:,1:3);
end

function R = rotxd(deg)
a = deg2rad(deg);
R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end

function Tstl = stlFrameFromPCA(V)
c = mean(V,1);
[~,~,Vh] = svd(V - c, 'econ');
R = Vh;
if det(R) < 0, R(:,3) = -R(:,3); end
Tstl = eye(4);
Tstl(1:3,1:3) = R;
Tstl(1:3,4)   = c(:);
end

function bb = bboxFromV(V)
bb = [min(V(:,1)) max(V(:,1)); ...
      min(V(:,2)) max(V(:,2)); ...
      min(V(:,3)) max(V(:,3))];
end

function N = faceNormalsFromFV(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
N = cross(v2 - v1, v3 - v1, 2);
nrm = sqrt(sum(N.^2,2));
nrm(nrm==0) = 1;
N = N ./ nrm;
end
