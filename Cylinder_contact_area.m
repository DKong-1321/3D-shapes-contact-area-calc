% USER OPTIONS
opts.trialIndices         = 1;           % [] = all trials
opts.animMode             = 'global';    % 'global' or 'tibiaFrame'
opts.contactTolerance_mm  = 0.10;

opts.viewAzEl   = [30 20];
opts.tibiaAlpha = 0.25;
opts.femurAlpha = 0.95;

% Colours
opts.tibiaColor   = [0.0 0.2 0.7];
opts.femurBaseRGB = [0.85 0.85 0.85];
opts.contactRGB   = [1.00 0.00 0.00];

% LOAD
projectRoot = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(projectRoot,'src')));
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','lib')));
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','external')));

S = load(fullfile(projectRoot,'data','matlab 1.mat'));
kinematics = S.kinematics;

femur_mesh = stlread(fullfile(projectRoot,'model','cylinder','top_matlab1.stl'));
tibia_mesh = stlread(fullfile(projectRoot,'model','cylinder','bottom_matlab1.stl'));

V_femur = femur_mesh.Points;
F_femur = femur_mesh.ConnectivityList;

V_tibia = tibia_mesh.Points;
F_tibia = tibia_mesh.ConnectivityList;

% TRIAL SELECTION
if isempty(opts.trialIndices)
    trialList = 1:numel(kinematics.trajectories);
else
    trialList = opts.trialIndices;
end

% MAIN LOOP
for TRIAL_INDEX = trialList

    tr = kinematics.trajectories(TRIAL_INDEX);
    cond = tr.LoadingCondition;
    if iscell(cond); cond = cond{1}; end
    cond = char(string(cond));

    gTfi = tr.Transform.gTfi;
    gTti = tr.Transform.gTti;
    fTt  = tr.Transform.fTt;
    nFrames = size(fTt,3);

    fprintf('\nRunning trial %d: %s\n', TRIAL_INDEX, cond);

    % CONTACT COMPUTE (ALWAYS IN TIBIA FRAME)
    tol = opts.contactTolerance_mm;

    tTf_all = pageinv(fTt);                 % ^tT_f  (femur into tibia frame)
    master  = buildBodyStruct(F_tibia, V_tibia);

    contactArea    = zeros(nFrames,1);
    contactMaskAll = false(size(F_femur,1), nFrames);

    for k = 1:nFrames
        Vf_k = applyRigid(tTf_all(:,:,k), V_femur);   % femur expressed in tibia frame
        slave = buildBodyStruct(F_femur, Vf_k);

        res = computeContactArea_STS_hybrid(master, slave, tol);

        contactArea(k)      = res.contactArea;
        contactMaskAll(:,k) = res.contactMask;
    end

    % PREP TRANSFORMS FOR RENDER
    V_femur_h = [V_femur'; ones(1,size(V_femur,1))];
    V_tibia_h = [V_tibia'; ones(1,size(V_tibia,1))];

    switch lower(opts.animMode)
        case 'global'
            Vf_all = pagemtimes(gTfi, V_femur_h);     % femur in global
            Vt_all = pagemtimes(gTti, V_tibia_h);     % tibia in global
        case 'tibiaframe'
            Vf_all = pagemtimes(tTf_all, V_femur_h);  % femur in tibia frame
            Vt_all = V_tibia;                         % tibia fixed
        otherwise
            error('opts.animMode must be ''global'' or ''tibiaFrame''.');
    end

    % ANIMATION (CONTACT OVERLAY)
    figure('Name',[cond ' — Contact Animation'],'Color','w');
    hold on; grid on; axis equal;
    view(opts.viewAzEl(1), opts.viewAzEl(2));
    camlight; lighting gouraud;

    % Tibia
    hTibia = patch('Faces',F_tibia,'Vertices',Vt_all(1:3,:,1)', ...
        'FaceColor',opts.tibiaColor, 'EdgeColor','none', 'FaceAlpha',opts.tibiaAlpha);

    % Femur (flat face colours => FaceVertexCData must be Nfaces x 1)
    hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)', ...
        'FaceColor','flat', 'EdgeColor','none', 'FaceAlpha',opts.femurAlpha);

    % 1 = base grey, 2 = contact red
    colormap([opts.femurBaseRGB; opts.contactRGB]);
    caxis([1 2]);
    C = ones(size(F_femur,1),1);

    for k = 1:nFrames
        % Update positions
        set(hFemur,'Vertices',Vf_all(1:3,:,k)');

        if strcmpi(opts.animMode,'global')
            set(hTibia,'Vertices',Vt_all(1:3,:,k)');
        end

        % Contact colouring
        C(:) = 1;
        C(contactMaskAll(:,k)) = 2;
        set(hFemur,'FaceVertexCData',C);

        title(sprintf('%s | Frame %d/%d | Contact Area = %.3f mm^2', ...
            cond, k, nFrames, contactArea(k)));

        drawnow;
    end

    hold off;

    % CONTACT AREA VS FRAME
    figure('Name',[cond ' — Contact Area vs Frame'],'Color','w');
    plot(contactArea,'LineWidth',1.5);
    xlabel('Frame');
    ylabel('Contact area (mm^2)');
    title(cond);
    grid on;

end

fprintf('\nDone.\n');

% HELPERS

function Vout = applyRigid(T,V)
Vh = [V ones(size(V,1),1)];
Vt = (T*Vh')';
Vout = Vt(:,1:3);
end

function body = buildBodyStruct(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
body.F = F;
body.V = V;
body.centroids = (v1+v2+v3)/3;
body.areas = 0.5*vecnorm(cross(v2-v1,v3-v1,2),2,2);
body.faceNormals = faceNormalsFromFV(F,V);
end

function N = faceNormalsFromFV(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
N = cross(v2 - v1, v3 - v1, 2);
mag = sqrt(sum(N.^2,2));
mag(mag==0) = 1;
N = N ./ mag;
end

% CONTACT AREA: HYBRID STS

function results = computeContactArea_STS_hybrid(master, slave, tol, opts)
if nargin < 3 || isempty(tol), tol = 0.2; end
if nargin < 4, opts = struct(); end 

master = ensureBodyFields_fast(master);
slave  = ensureBodyFields_fast(slave);

if ~isfield(master,'faceNormals') || isempty(master.faceNormals)
    master.faceNormals = faceNormalsFromFV(master.F, master.V);
end
if ~isfield(master,'kdtreeFaces') || isempty(master.kdtreeFaces)
    master.kdtreeFaces = KDTreeSearcher(master.centroids);
end

[Apen, Atol, maskPen, maskTol] = computeHybridArea(master, slave.F, slave.V, tol);

results = struct();
results.contactArea     = Atol;     % tolerance-band clipped area
results.penetrationArea = Apen;     % penetration-only clipped area
results.contactMask     = maskTol;  % per-slave-face mask
results.penMask         = maskPen;
results.tolMask         = maskTol;
end

function [Apen, Atol, maskPen, maskTol] = computeHybridArea(master, Fslave, Vslave, tol)
P1 = Vslave(Fslave(:,1),:);
P2 = Vslave(Fslave(:,2),:);
P3 = Vslave(Fslave(:,3),:);
C  = (P1 + P2 + P3) / 3;

idx = knnsearch(master.kdtreeFaces, C, 'K', 1);

Fm  = master.F(idx,:);
V1m = master.V(Fm(:,1),:);
Nm  = master.faceNormals(idx,:);

d1 = dot(P1 - V1m, Nm, 2);
d2 = dot(P2 - V1m, Nm, 2);
d3 = dot(P3 - V1m, Nm, 2);

Aclip0   = zeros(size(Fslave,1),1);
AclipTol = zeros(size(Fslave,1),1);

for i = 1:size(Fslave,1)
    Aclip0(i)   = clippedAreaTriangleHalfSpace(P1(i,:),P2(i,:),P3(i,:), d1(i),d2(i),d3(i));
    AclipTol(i) = clippedAreaTriangleHalfSpace(P1(i,:),P2(i,:),P3(i,:), d1(i)-tol, d2(i)-tol, d3(i)-tol);
end

Apen    = sum(Aclip0);
Atol    = sum(AclipTol);
maskPen = Aclip0   > 0;
maskTol = AclipTol > 0;
end

function A = clippedAreaTriangleHalfSpace(p1,p2,p3, d1,d2,d3)
pts = [p1; p2; p3];
ds  = [d1; d2; d3];

inside = ds <= 0;
if all(~inside), A = 0; return; end
if all(inside),  A = triArea3D(p1,p2,p3); return; end

poly = [];
for k = 1:3
    k2 = mod(k,3) + 1;
    Pk  = pts(k,:);  dk  = ds(k);
    Pk2 = pts(k2,:); dk2 = ds(k2);

    if dk <= 0
        poly = [poly; Pk]; 
    end

    if (dk <= 0 && dk2 > 0) || (dk > 0 && dk2 <= 0)
        t  = dk / (dk - dk2);
        Pi = Pk + t*(Pk2 - Pk);
        poly = [poly; Pi]; 
    end
end

if size(poly,1) < 3, A = 0; return; end

A = 0;
p0 = poly(1,:);
for j = 2:(size(poly,1)-1)
    A = A + triArea3D(p0, poly(j,:), poly(j+1,:));
end
end

function A = triArea3D(a,b,c)
A = 0.5 * norm(cross(b-a, c-a));
end

function body = ensureBodyFields_fast(body)
if ~isfield(body,'centroids')
    if isfield(body,'triCentroid')
        body.centroids = body.triCentroid;
    else
        error('Body missing centroids/triCentroid.');
    end
end
if ~isfield(body,'areas')
    if isfield(body,'triArea')
        body.areas = body.triArea;
    else
        error('Body missing areas/triArea.');
    end
end
if ~isfield(body,'bbox') || isempty(body.bbox)
    V = body.V;
    body.bbox = [min(V(:,1)) max(V(:,1)); min(V(:,2)) max(V(:,2)); min(V(:,3)) max(V(:,3))];
end
end