%% cube_cylinder_contact_animation_4x4.m
% Cube–cylinder contact animation
% Motion entirely via 4x4 homogeneous transforms
% Contact patch shown correctly

clear; clc; close all;

%% ---------------- PATHS ----------------
modelDir = fullfile(pwd,'model');

cubeSTL = fullfile(modelDir,'50mm cube 4 iterations.stl');
cylSTL  = fullfile(modelDir,'cylinder 30mm diameter 4 it.stl');

assert(exist(cubeSTL,'file')==2, 'Missing cube STL');
assert(exist(cylSTL,'file')==2,  'Missing cylinder STL');

%% ---------------- OPTIONS ----------------
tol         = 0.3;     % mm
penetration = 0.8;     % mm

nFrames    = 150;
framePause = 0.02;

slideAmp   = 12;       % mm
yawAmpDeg  = 25;       % deg

%% ---------------- LOAD STLS ----------------
Sc = stlread(cubeSTL);
Sy = stlread(cylSTL);

Vc  = Sc.Points;  Fc = Sc.ConnectivityList;
Vy0 = Sy.Points;  Fy = Sy.ConnectivityList;

%% ---------------- BUILD MASTER (CUBE) ----------------
master = buildBodyStruct(Fc, Vc);

%% ---------------- INITIAL PLACEMENT ----------------
cubeTopZ = max(Vc(:,3));
cCube    = mean(Vc,1);
cCyl0    = mean(Vy0,1);

Txy = eye(4);
Txy(1:3,4) = (cCube - cCyl0)';

Vy_xy = applyTform4x4(Txy, Vy0);
cylBottomZ = min(Vy_xy(:,3));

Tz = eye(4);
Tz(3,4) = cubeTopZ - cylBottomZ - penetration;

Tbase = Tz * Txy;

%% ---------------- FIGURE ----------------
figure('Color','w'); hold on; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');

% Cube (master)
patch('Faces',Fc,'Vertices',Vc, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'FaceAlpha',0.25, ...
      'EdgeColor','none');

% Cylinder (slave)
hCyl = patch('Faces',Fy,'Vertices',Vy0, ...
             'FaceColor','flat', ...
             'FaceAlpha',0.5, ...
             'EdgeColor','none');

% Scalar colormap:
% 1 = non-contact, 2 = contact
baseColor    = [0.92 0.92 0.92];
contactColor = [0.0  1.0  0.2];

colormap([baseColor;
          contactColor]);
clim([1 2]);

camlight headlight;
lighting gouraud;

%% ---------------- ANIMATION ----------------
for k = 1:nFrames

    t = (k-1)/(nFrames-1);

    % Translation (X)
    Tx = eye(4);
    Tx(1,4) = slideAmp * cos(2*pi*t);

    % Rotation about Z (yaw)
    yaw = deg2rad(yawAmpDeg * sin(2*pi*t));
    Rz = eye(4);
    Rz(1:3,1:3) = [cos(yaw) -sin(yaw) 0;
                   sin(yaw)  cos(yaw) 0;
                   0         0        1];

    % Rotate about cylinder centroid
    Vy_base = applyTform4x4(Tbase, Vy0);
    cBase   = mean(Vy_base,1);

    Tto0 = eye(4);  Tto0(1:3,4) = -cBase(:);
    Tbk  = eye(4);  Tbk(1:3,4)  =  cBase(:);

    % Full 4x4 transform
    T = Tx * Tbk * Rz * Tto0 * Tbase;

    % Apply transform
    Vy = applyTform4x4(T, Vy0);

    slave = buildBodyStruct(Fy, Vy);

    % Contact
    res = computeContactArea_STS(master, slave, tol);

    % ---- PER-FACE SCALAR COLOURING ----
    C = ones(size(Fy,1),1);      % non-contact
    if ~isempty(res.contactMask)
        C(res.contactMask) = 2;  % contact
    end

    set(hCyl,'Vertices',Vy,'CData',C);

    title(sprintf('Frame %d / %d | Contact area = %.3f mm^2 | Contact tris = %d', ...
        k, nFrames, res.contactArea, nnz(res.contactMask)));

    drawnow;
    pause(framePause);
end

%% ---------------- HELPERS ----------------
function body = buildBodyStruct(F,V)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);
    body.F = F;
    body.V = V;
    body.centroids = (v1+v2+v3)/3;
    body.areas = 0.5 * vecnorm(cross(v2-v1,v3-v1,2),2,2);
end

function Vout = applyTform4x4(T,V)
    Vh = [V ones(size(V,1),1)];
    Vt = (T * Vh')';
    Vout = Vt(:,1:3);
end

