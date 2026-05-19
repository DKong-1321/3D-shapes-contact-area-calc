% run_cylinder_visualisation_contactVis.m
% Place in root of 3D-shapes-contact-area-calc and run from that folder.
% Requires: data/matlab 1.mat, model/cylinder/top_matlab1.stl, bottom_matlab1.stl
%
% Adds:
%   - Trial selection
%   - Global or tibia-frame STL animation
%   - Contact area computation (in tibia frame)
%   - Contact mask visualisation on femur faces (grey vs green)
%   - Optional plots

clear; clc; close all;

%% ========================== USER OPTIONS ==========================
opts = struct();

% --- What to run ---
opts.doDigitisationFigure   = false;   % Figure 1
opts.doAlignmentPreview     = false;   % Figure 2 (frame 1 sanity check)
opts.doBodyAxesAnimation    = false;   % axes animation
opts.doKinematicsPlots      = false;   % 6 metrics (3x2)
opts.doSTLAnimation         = true;    % STL animation window

% --- Contact ---
opts.doContactCompute       = true;    % compute contact area vs frame
opts.doContactVisualisation = true;    % colour contact tris on femur (requires doContactCompute)
opts.contactTolerance_mm    = 0.20;    % mm
opts.showContactText        = true;    % print area on title
opts.plotContactVsFrame     = true;
opts.plotContactVsFlexion   = false;   % only meaningful if doKinematicsPlots = true or you compute flexion anyway

% --- Trials to run ---
%   []    -> all trials
%   5     -> just 5th run
%   [2 5] -> run 2 and 5
opts.trialIndices = 5;

% --- Animation mode ---
%   'tibiaFrame' -> tibia fixed, femur moved by tTf = inv(fTt)
%   'global'     -> BOTH move in global using gTfi/gTti (tibia rotation visible)
opts.animMode = 'global';

% --- Video saving ---
opts.saveVideo  = false;    % default OFF
opts.frameRate  = 30;

% --- Visual settings ---
opts.viewAzEl   = [30 20];
opts.faceAlpha  = 0.85;
opts.pausePerFrame = 0.0;   % set >0 to slow down (e.g. 0.01)

% Colours
opts.tibiaColor      = [0.55 0.70 0.85];
opts.femurBaseColor  = [0.75 0.75 0.75];   % non-contact
opts.contactColor    = [0.00 1.00 0.20];   % contact (bright green)

% Alignment preview uses this trial
opts.previewTrial = 1;

%% ===================== PROJECT SETUP / PATHS ======================
projectRoot = fileparts(mfilename('fullpath'));

% Your contact code lives under src/
addpath(genpath(fullfile(projectRoot,'src')));

% OpticalTracking classes
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','lib')));
addpath(genpath(fullfile(projectRoot,'OpticalTracking-dev-3','external')));

matFile = fullfile(projectRoot,'data','matlab 1.mat');
assert(exist(matFile,'file')==2,'Cannot find data/matlab 1.mat');
S = load(matFile);

digitisation = S.digitisation;
kinematics   = S.kinematics;
config       = digitisation(1).config;

nTrials = numel(kinematics.trajectories);
fprintf('Loaded trajectories: %d\n', nTrials);

% Trial list
if isempty(opts.trialIndices)
    trialList = 1:nTrials;
else
    trialList = opts.trialIndices(:)';
end

% STL paths
femurStlPath = fullfile(projectRoot,'model','cylinder','top_matlab1.stl');
tibiaStlPath = fullfile(projectRoot,'model','cylinder','bottom_matlab1.stl');
assert(exist(femurStlPath,'file')==2,'Missing: %s',femurStlPath);
assert(exist(tibiaStlPath,'file')==2,'Missing: %s',tibiaStlPath);

% Load STLs (shared across all trials)
femur_mesh  = stlread(femurStlPath);
tibia_mesh  = stlread(tibiaStlPath);

V_femur = femur_mesh.Points;
F_femur = femur_mesh.ConnectivityList;

V_tibia = tibia_mesh.Points;
F_tibia = tibia_mesh.ConnectivityList;

fprintf('STL: femur faces=%d verts=%d | tibia faces=%d verts=%d\n', ...
    size(F_femur,1), size(V_femur,1), size(F_tibia,1), size(V_tibia,1));

%% ====================== FIGURE 1: Digitisation =====================
if opts.doDigitisationFigure
    digitisation.visualise();
end

%% ================== FIGURE 2: Alignment Preview ====================
if opts.doAlignmentPreview
    pTrial = min(max(opts.previewTrial,1), nTrials);
    tr0 = kinematics.trajectories(pTrial);

    gTfi0 = tr0.Transform.gTfi;
    gTti0 = tr0.Transform.gTti;
    fTt0  = tr0.Transform.fTt;

    figure('Name','Alignment preview (frame 1)','Color','w');
    hold on; grid on; axis equal; view(opts.viewAzEl(1),opts.viewAzEl(2));
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

    switch lower(opts.animMode)
        case 'tibiaframe'
            tTf1 = inv(fTt0(:,:,1)); % femur in tibia frame
            patch('Faces',F_tibia,'Vertices',V_tibia,'EdgeColor','none', ...
                'FaceColor',opts.tibiaColor,'FaceAlpha',0.30);
            patch('Faces',F_femur,'Vertices',applyRigid(tTf1,V_femur),'EdgeColor','none', ...
                'FaceColor',[0.9 0.7 0.6],'FaceAlpha',0.30);
            title(sprintf('Preview: tibia frame (trial %d, frame 1)', pTrial));

        case 'global'
            Vt_g = applyRigid(gTti0(:,:,1), V_tibia);
            Vf_g = applyRigid(gTfi0(:,:,1), V_femur);
            patch('Faces',F_tibia,'Vertices',Vt_g,'EdgeColor','none', ...
                'FaceColor',opts.tibiaColor,'FaceAlpha',0.30);
            patch('Faces',F_femur,'Vertices',Vf_g,'EdgeColor','none', ...
                'FaceColor',[0.9 0.7 0.6],'FaceAlpha',0.30);
            title(sprintf('Preview: global frame (trial %d, frame 1)', pTrial));

        otherwise
            error('opts.animMode must be ''tibiaFrame'' or ''global''.');
    end

    camlight; lighting gouraud;
    hold off;
end

%% ========================= LOOP OVER TRIALS =========================
for TRIAL_INDEX = trialList

    tr   = kinematics.trajectories(TRIAL_INDEX);
    cond = tr.LoadingCondition;

    gTfi = tr.Transform.gTfi;
    gTti = tr.Transform.gTti;
    fTt  = tr.Transform.fTt;

    nFrames = size(fTt,3);

    fprintf('\n── Trial %d / %d : %s | frames=%d ──\n', TRIAL_INDEX, nTrials, cond, nFrames);

    %% ---------- Body axes animation ----------
    if opts.doBodyAxesAnimation
        femur_origins = squeeze(gTfi(1:3,4,:));
        tibia_origins = squeeze(gTti(1:3,4,:));

        SCALE = 30; STEP = 3;
        figAxes = figure('Name', sprintf('Body frame axes — %s', cond), 'Color','w');

        for k = 1:STEP:nFrames
            clf(figAxes);
            hold on; grid on; axis equal; view(opts.viewAzEl(1),opts.viewAzEl(2));
            xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
            title(sprintf('%s — frame %d / %d', cond, k, nFrames));

            Rf = gTfi(1:3,1:3,k); of = gTfi(1:3,4,k);
            Rt = gTti(1:3,1:3,k); ot = gTti(1:3,4,k);

            for ax = 1:3
                quiver3(of(1),of(2),of(3), SCALE*Rf(1,ax),SCALE*Rf(2,ax),SCALE*Rf(3,ax), 0,'r','LineWidth',1.5);
                quiver3(ot(1),ot(2),ot(3), SCALE*Rt(1,ax),SCALE*Rt(2,ax),SCALE*Rt(3,ax), 0,'b','LineWidth',1.5);
            end

            plot3(femur_origins(1,1:k),femur_origins(2,1:k),femur_origins(3,1:k),'r:','LineWidth',0.5);
            plot3(tibia_origins(1,1:k),tibia_origins(2,1:k),tibia_origins(3,1:k),'b:','LineWidth',0.5);

            drawnow;
            pause(0.03);
            hold off;
        end
    end

    %% ---------- Kinematics (all six metrics) ----------
    tf = [];
    if opts.doKinematicsPlots || opts.plotContactVsFlexion
        bones.femur.i      = squeeze(gTfi(1:3,1,:));
        bones.femur.j      = squeeze(gTfi(1:3,2,:));
        bones.femur.k      = squeeze(gTfi(1:3,3,:));
        bones.femur.origin = squeeze(gTfi(1:3,4,:));
        bones.tibia.i      = squeeze(gTti(1:3,1,:));
        bones.tibia.j      = squeeze(gTti(1:3,2,:));
        bones.tibia.k      = squeeze(gTti(1:3,3,:));
        bones.tibia.origin = squeeze(gTti(1:3,4,:));

        [tf, ~] = Knee.grood_and_suntay(bones.femur, bones.tibia, [], fTt, [], config.is_right_knee);
    end

    if opts.doKinematicsPlots
        metrics = {'flexion',  'Flexion (deg)'; ...
                   'varus',    'Varus (deg)'; ...
                   'external', 'External rotation (deg)'; ...
                   'lateral',  'Lateral (mm)'; ...
                   'anterior', 'Anterior (mm)'; ...
                   'superior', 'Superior (mm)'};

        figure('Name', sprintf('Kinematics — %s', cond), 'Color','w');
        for m = 1:6
            subplot(3,2,m);
            plot(1:nFrames, tf.(metrics{m,1}), 'k', 'LineWidth', 1.2);
            xlabel('Frame'); ylabel(metrics{m,2});
            title(metrics{m,2});
            grid on;
        end
        sgtitle(sprintf('Tibiofemoral kinematics — %s', cond));
    end

    %% ---------- Contact compute (always in tibia frame) ----------
    contactArea = [];
    contactMaskAll = [];

    if opts.doContactCompute || opts.doContactVisualisation

        tol = opts.contactTolerance_mm;

        % Femur in tibia frame: tTf(:,:,k) = inv(fTt(:,:,k))
        tTf = pageinv(fTt);

        master = buildBodyStruct(F_tibia, V_tibia);

        contactArea = zeros(nFrames,1);

        if opts.doContactVisualisation
            contactMaskAll = false(size(F_femur,1), nFrames); % per-face per-frame
        end

        for k = 1:nFrames
            Vf_k = applyRigid(tTf(:,:,k), V_femur);
            slave = buildBodyStruct(F_femur, Vf_k);

            res = computeContactArea_STS_hybrid(master, slave, tol);
            contactArea(k) = res.contactArea;

            if opts.doContactVisualisation
                if isfield(res,'contactMask') && numel(res.contactMask)==size(F_femur,1)
                    contactMaskAll(:,k) = logical(res.contactMask(:));
                else
                    % If your hybrid function is configured not to return a mask, you'll just see no highlighting.
                    contactMaskAll(:,k) = false(size(F_femur,1),1);
                end
            end
        end

        fprintf('Contact: mean=%.4f | max=%.4f (mm^2)\n', mean(contactArea), max(contactArea));

        if opts.plotContactVsFrame
            figure('Name', sprintf('Contact area vs frame — %s', cond), 'Color','w');
            plot(1:nFrames, contactArea, 'LineWidth', 1.5);
            xlabel('Frame'); ylabel('Contact area (mm^2)');
            title(sprintf('Contact area vs frame — %s', cond));
            grid on;
        end

        if opts.plotContactVsFlexion && ~isempty(tf)
            figure('Name', sprintf('Contact area vs flexion — %s', cond), 'Color','w');
            plot(tf.flexion, contactArea, 'LineWidth', 1.5);
            xlabel('Flexion (deg)'); ylabel('Contact area (mm^2)');
            title(sprintf('Contact area vs flexion — %s', cond));
            grid on;
        end
    end

    %% ---------- STL animation (with optional contact colouring) ----------
    if opts.doSTLAnimation

        figStl = figure('Name', sprintf('STL — %s', cond), 'Color', 'w');
        clf(figStl);
        hold on; grid on; axis equal; view(opts.viewAzEl(1), opts.viewAzEl(2));
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

        % Precompute global or tibia-frame vertex arrays for speed
        V_femur_h = [V_femur'; ones(1,size(V_femur,1))];
        V_tibia_h = [V_tibia'; ones(1,size(V_tibia,1))];

        switch lower(opts.animMode)

            case 'tibiaframe'
                % Femur moves in tibia frame
                tTf = pageinv(fTt);
                Vf_all = pagemtimes(tTf, V_femur_h);      % [4 x Nv_f x N]
                % Tibia fixed
                Vt_fixed = V_tibia;

                % Axis ranges
                all_femur = reshape(Vf_all(1:3,:,:), 3, [])';
                all_verts = [all_femur; Vt_fixed];
                x_range = [min(all_verts(:,1)) max(all_verts(:,1))];
                y_range = [min(all_verts(:,2)) max(all_verts(:,2))];
                z_range = [min(all_verts(:,3)) max(all_verts(:,3))];

                hTibia = patch('Faces',F_tibia,'Vertices',Vt_fixed, ...
                    'FaceColor',opts.tibiaColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

                % Femur patch (flat per-face colour)
                hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)', ...
                    'FaceColor','flat','EdgeColor','none','FaceAlpha',opts.faceAlpha, ...
                    'CDataMapping','direct');

            case 'global'
                % Both move in global
                Vf_all = pagemtimes(gTfi, V_femur_h);
                Vt_all = pagemtimes(gTti, V_tibia_h);

                all_femur = reshape(Vf_all(1:3,:,:), 3, [])';
                all_tibia = reshape(Vt_all(1:3,:,:), 3, [])';
                all_verts = [all_femur; all_tibia];

                x_range = [min(all_verts(:,1)) max(all_verts(:,1))];
                y_range = [min(all_verts(:,2)) max(all_verts(:,2))];
                z_range = [min(all_verts(:,3)) max(all_verts(:,3))];

                hTibia = patch('Faces',F_tibia,'Vertices',Vt_all(1:3,:,1)', ...
                    'FaceColor',opts.tibiaColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

                hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)', ...
                    'FaceColor','flat','EdgeColor','none','FaceAlpha',opts.faceAlpha, ...
                    'CDataMapping','direct');

            otherwise
                error('opts.animMode must be ''tibiaFrame'' or ''global''.');
        end

        xlim(x_range); ylim(y_range); zlim(z_range);
        camlight; lighting gouraud;
        legend([hFemur hTibia], {'Femur','Tibia'});

        % Two-colour map: 1=base (grey), 2=contact (green)
        colormap([opts.femurBaseColor; opts.contactColor]);
        caxis([1 2]);

        % Optional video writer
        if opts.saveVideo
            safeName = regexprep(cond,'[^a-zA-Z0-9_-]','_');
            videoFile = fullfile(projectRoot, sprintf('cylinder_%s.mp4', safeName));
            vw = VideoWriter(videoFile,'MPEG-4');
            vw.FrameRate = opts.frameRate;
            open(vw);
        end

        % Initial colour data
        C = ones(size(F_femur,1),1);

        for k = 1:nFrames

            % Update vertices
            set(hFemur,'Vertices',Vf_all(1:3,:,k)');

            if strcmpi(opts.animMode,'global')
                set(hTibia,'Vertices',Vt_all(1:3,:,k)');
            end

            % Update contact colouring
            if opts.doContactVisualisation && ~isempty(contactMaskAll)
                mask = contactMaskAll(:,k);
                C(:) = 1;
                C(mask) = 2;
                set(hFemur,'FaceVertexCData',C);
            else
                set(hFemur,'FaceVertexCData',C); % all grey
            end

            % Title (optionally include area)
            if opts.showContactText && ~isempty(contactArea)
                title(sprintf('%s — frame %d/%d | A=%.4f mm^2', cond, k, nFrames, contactArea(k)));
            else
                title(sprintf('%s — frame %d/%d', cond, k, nFrames));
            end

            drawnow;

            if opts.saveVideo
                writeVideo(vw, getframe(figStl));
            end

            if opts.pausePerFrame > 0
                pause(opts.pausePerFrame);
            end
        end

        if opts.saveVideo
            close(vw);
            fprintf('Saved: %s\n', videoFile);
        end

        hold off;
    end
end

fprintf('\nAll selected trials done.\n');

%% ============================ HELPERS =============================
function Vout = applyRigid(T, V)
    Vh = [V ones(size(V,1),1)];
    Vt = (T * Vh')';
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
end