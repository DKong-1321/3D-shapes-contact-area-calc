%OPTIONS

opts.trialIndices = [];          % [] = all trials
opts.animMode     = 'global';   % 'global' or 'tibiaFrame'
opts.contactTolerance_mm = 0.50;

opts.viewAzEl = [30 20];

% Visual
opts.tibiaAlpha = 0.25;
opts.femurAlpha = 0.95;

% Mesh refinement options 
opts.refine.enable      = true;

% Choose ONE:
%   'neutralProximity' : frame-1 proximity ROI (distance threshold)
%   'contactSeed'      : run contact once (frame-1), refine union contact region
opts.refine.method      = 'contactSeed';

opts.refine.levels      = 2;      % how many refinement passes (1–3 is typical)
opts.refine.expandRings = 2;      % expand ROI by adjacency rings each level

% For neutralProximity
opts.refine.proxDist_mm = 1.0;    % ROI distance threshold in tibia frame

% For contactSeed
opts.refine.seedTol_mm  = 0.50;   % seed contact tolerance (bigger than final tol)

% Save animation to MP4 
opts.saveVideo     = true;   % set false if you don't want mp4
opts.videoFPS      = 30;     % frames per second in output video
opts.videoQuality  = 95;     % 0–100
% opts.videoFolder will be set after projectRoot is known

% LOAD

projectRoot = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(projectRoot,'src')));
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','lib')));
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','external')));

opts.videoFolder = fullfile(projectRoot,'results'); % default save folder

S = load(fullfile(projectRoot,'data','matlab 1.mat'));
kinematics = S.kinematics;

femur_mesh = stlread(fullfile(projectRoot,'model','cylinder','top_matlab1.stl'));
tibia_mesh = stlread(fullfile(projectRoot,'model','cylinder','bottom_matlab1.stl'));

V_femur0 = femur_mesh.Points;
F_femur0 = femur_mesh.ConnectivityList;
V_tibia0 = tibia_mesh.Points;
F_tibia0 = tibia_mesh.ConnectivityList;

%% ================= TRIAL SELECTION =================

if isempty(opts.trialIndices)
    trialList = 1:numel(kinematics.trajectories);
else
    trialList = opts.trialIndices;
end

%% ================= MAIN LOOP =================

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

    % Reset meshes each trial (so refinement doesn't accumulate across trials)
    F_femur = F_femur0; V_femur = V_femur0;
    F_tibia = F_tibia0; V_tibia = V_tibia0;

    % OPTIONAL: PRE-REFINE MESHES
    % Determine an ROI using frame 1 in tibia frame.
    % Then refine BOTH tibia and femur around the ROI.

    tTf_all = pageinv(fTt); % femur into tibia frame transforms
    tTf1 = tTf_all(:,:,1);

    if opts.refine.enable
        fprintf('Refinement: %s | levels=%d | expandRings=%d\n', ...
            opts.refine.method, opts.refine.levels, opts.refine.expandRings);

        for lvl = 1:opts.refine.levels

            % Build pose for ROI selection (frame 1 in tibia frame)
            Vf1 = applyRigid(tTf1, V_femur);

            master0 = buildBodyStruct(F_tibia, V_tibia);
            slave0  = buildBodyStruct(F_femur, Vf1);

            switch lower(opts.refine.method)

                case 'neutralproximity'
                    % Mark femur faces whose centroid is within proxDist to nearest tibia face centroid
                    idxM = knnsearch(master0.centroids, slave0.centroids);
                    d = vecnorm(slave0.centroids - master0.centroids(idxM,:), 2, 2);

                    maskF = d < opts.refine.proxDist_mm;

                    % Also mark corresponding tibia faces (nearest ones)
                    maskT = false(size(F_tibia,1),1);
                    maskT(idxM(maskF)) = true;

                case 'contactseed'
                    % Run contact once with a "fat" tolerance to seed the ROI
                    resSeed = computeContactArea_STS_hybrid(master0, slave0, opts.refine.seedTol_mm);

                    maskF = resSeed.contactMask;

                    % Mark corresponding tibia faces via nearest-centroid mapping
                    idxM = knnsearch(master0.centroids, slave0.centroids);
                    maskT = false(size(F_tibia,1),1);
                    maskT(idxM(maskF)) = true;

                otherwise
                    error('Unknown refine.method: %s', opts.refine.method);
            end

            % Expand ROI by adjacency rings
            for r = 1:opts.refine.expandRings
                maskF = expandOneRing(F_femur, maskF);
                maskT = expandOneRing(F_tibia, maskT);
            end

            if ~any(maskF) && ~any(maskT)
                warning('No ROI faces found at refine level %d — stopping refinement.', lvl);
                break
            end

            % Refine BOTH meshes around ROI
            [F_femur, V_femur] = refineTrianglesMidpoint_cached_simple(F_femur, V_femur, maskF);
            [F_tibia, V_tibia] = refineTrianglesMidpoint_cached_simple(F_tibia, V_tibia, maskT);

            fprintf('  refine lvl %d done | femur faces=%d | tibia faces=%d\n', ...
                lvl, size(F_femur,1), size(F_tibia,1));
        end
    end

    % CONTACT COMPUTE (FINAL MESH)

    contactArea    = zeros(nFrames,1);
    contactMaskAll = false(size(F_femur,1), nFrames);

    % Compute contact in tibia frame (tibia fixed; femur transformed each frame)
    master = buildBodyStruct(F_tibia, V_tibia);

    for k = 1:nFrames
        Vf_k = applyRigid(tTf_all(:,:,k), V_femur);
        slave = buildBodyStruct(F_femur, Vf_k);

        res = computeContactArea_STS_hybrid(master, slave, opts.contactTolerance_mm);

        contactArea(k) = res.contactArea;

        % res.contactMask is per FEMUR face (works because slave.F is femur F)
        contactMaskAll(:,k) = res.contactMask;
    end

    % PREP TRANSFORMS FOR RENDER

    V_femur_h = [V_femur'; ones(1,size(V_femur,1))];
    V_tibia_h = [V_tibia'; ones(1,size(V_tibia,1))];

    switch lower(opts.animMode)

        case 'global'
            % Note: gTfi/gTti come from tracking; they still apply to the geometry.
            Vf_all = pagemtimes(gTfi, V_femur_h);
            Vt_all = pagemtimes(gTti, V_tibia_h);

        case 'tibiaframe'
            Vf_all = pagemtimes(tTf_all, V_femur_h);
            Vt_all = V_tibia; % fixed

        otherwise
            error('Unknown animMode: %s', opts.animMode);
    end

    % CONTACT ANIMATION

    figAnim = figure('Name',[cond ' — Contact Animation'],'Color','w');
    set(figAnim,'Renderer','opengl'); % better for patch + lighting capture
    hold on; grid on; axis equal;
    view(opts.viewAzEl);
    camlight; lighting gouraud;

    hTibia = patch('Faces',F_tibia,'Vertices',Vt_all(1:3,:,1)',...
        'FaceColor',[0 0.2 0.7],...
        'EdgeColor','none',...
        'FaceAlpha',opts.tibiaAlpha);

    hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)',...
        'FaceColor','flat',...
        'EdgeColor','none',...
        'FaceAlpha',opts.femurAlpha);

    % Grey + Red
    colormap([0.85 0.85 0.85;
              1.00 0.00 0.00]);
    caxis([1 2]);

    C = ones(size(F_femur,1),1);

    % Video writer setup (per trial) 
    if opts.saveVideo
        if ~exist(opts.videoFolder,'dir')
            mkdir(opts.videoFolder);
        end

        safeCond = regexprep(cond,'[^a-zA-Z0-9_-]','_');
        mp4Path  = fullfile(opts.videoFolder, sprintf('trial_%02d_%s_contact.mp4', TRIAL_INDEX, safeCond));

        vw = VideoWriter(mp4Path, 'MPEG-4');
        vw.FrameRate = opts.videoFPS;
        vw.Quality   = opts.videoQuality;
        open(vw);

        fprintf('Saving MP4 to: %s\n', mp4Path);
    end

    for k = 1:nFrames

        set(hFemur,'Vertices',Vf_all(1:3,:,k)');

        if strcmpi(opts.animMode,'global')
            set(hTibia,'Vertices',Vt_all(1:3,:,k)');
        end

        C(:) = 1;
        C(contactMaskAll(:,k)) = 2;

        set(hFemur,'FaceVertexCData',C);

        title(sprintf('%s | Frame %d | Contact Area = %.3f mm^2', ...
            cond,k,contactArea(k)));

        drawnow;

        if opts.saveVideo
            fr = getframe(figAnim);
            writeVideo(vw, fr);
        end
    end

    if opts.saveVideo
        close(vw);
        fprintf('Video written.\n');
    end

    hold off;

    % CONTACT AREA VS FRAME

    figure('Name',[cond ' — Contact Area']);
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
N = cross(v2-v1,v3-v1,2);
mag = sqrt(sum(N.^2,2));
mag(mag==0)=1;
N = N./mag;
end

function mask2 = expandOneRing(F, mask)
mask2 = mask;

nV = max(F(:));
v2f = cell(nV,1);

for fi = 1:size(F,1)
    for j = 1:3
        v = F(fi,j);
        v2f{v}(end+1) = fi;
    end
end

for fi = find(mask)'
    verts = F(fi,:);
    for v = verts
        mask2(v2f{v}) = true;
    end
end
end

function [Fnew, Vnew] = refineTrianglesMidpoint_cached_simple(F, V, mask)
% Midpoint subdivision of selected triangles.
% Refines each marked triangle into 4 triangles, reusing mid-edge vertices.

Vnew = V;
nF = size(F,1);

% worst case: each refined face adds +3 faces (1 -> 4)
Fnew = zeros(nF + 3*nnz(mask), 3);

edgeMid = containers.Map('KeyType','char','ValueType','int32');

o = 0;
for fi = 1:nF
    a = F(fi,1); b = F(fi,2); c = F(fi,3);

    if mask(fi)
        ab = mid(a,b);
        bc = mid(b,c);
        ca = mid(c,a);

        o=o+1; Fnew(o,:)=[a  ab ca];
        o=o+1; Fnew(o,:)=[b  bc ab];
        o=o+1; Fnew(o,:)=[c  ca bc];
        o=o+1; Fnew(o,:)=[ab bc ca];
    else
        o=o+1; Fnew(o,:)=[a b c];
    end
end

Fnew = Fnew(1:o,:);

    function m = mid(i,j)
        key = sprintf('%d_%d', min(i,j), max(i,j));
        if isKey(edgeMid,key)
            m = edgeMid(key);
        else
            Vnew(end+1,:) = 0.5*(Vnew(i,:) + Vnew(j,:));
            m = size(Vnew,1);
            edgeMid(key) = m;
        end
    end
end