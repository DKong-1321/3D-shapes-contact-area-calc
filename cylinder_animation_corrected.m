% run_cylinder_visualisation.m
% Place in root of 3D-shapes-contact-area-calc and run from that folder.
% Requires: data/matlab 1.mat, model/cylinder/top_matlab1.stl, bottom_matlab1.stl

clear; clc; close all;

%% ========================== USER OPTIONS ==========================
opts = struct();

% What to run
opts.doDigitisationFigure   = true;    % Figure 1
opts.doAlignmentPreview     = true;    % Figure 2 (frame 1 check)
opts.doBodyAxesAnimation    = false;   % body frame axes anim
opts.doKinematicsPlots      = false;   % 6 metrics plot (3x2)
opts.doSTLAnimation         = true;    % STL animation figure
opts.saveVideo              = false;   % << set true only if you want mp4

% Which trials (empty = all)
opts.trialIndices = [];   % e.g. [1 3 5]

% Animation mode:
%   'tibiaFrame' -> tibia fixed, femur moves by tTf = inv(fTt)
%   'global'     -> both move in global using gTfi and gTti  (THIS is what you want)
opts.animMode = 'global';

% Visual settings
opts.viewAzEl    = [30 20];
opts.pauseAxes   = 0.03;
opts.stepAxes    = 3;
opts.frameRate   = 30;

% STL appearance
opts.femurColor  = [0.85 0.70 0.55];
opts.tibiaColor  = [0.55 0.70 0.85];
opts.faceAlpha   = 0.85;

% Alignment preview 
opts.previewTrial = 1;

%% ===================== ERROR INJECTION CONFIG =====================
ERR = struct();
ERR.active = false;

% Landmark errors (not used in this simplified script; kept for compatibility)
ERR.swap_medial_lateral  = false;
ERR.flip_distal_z        = false;
ERR.offset_medial_mm     = 0;

% STL errors (apply to femur only, same as your earlier script)
ERR.stl_rotate_90_about_x = false;
ERR.stl_flip_z            = false;
ERR.stl_offset_z_mm       = 0;

% Transform errors
ERR.fTt_add_flexion_deg  = 0;
ERR.fTt_swap_tibia_femur = false;

%% ============================== SETUP =============================
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'lib')));
addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'external')));

required = {'Digitisation', 'JCS', 'Kinematics', 'Trajectory'};
for i = 1:numel(required)
    if exist(required{i}, 'class') == 8, fprintf(' OK: %s\n', required{i});
    else, fprintf(' MISSING: %s\n', required{i}); end
end

matFile = fullfile(projectRoot, 'data', 'matlab 1.mat');
assert(exist(matFile, 'file') == 2, 'Cannot find data/matlab 1.mat');

S = load(matFile);
digitisation = S.digitisation;
kinematics   = S.kinematics;
config       = digitisation(1).config;

nTrials = numel(kinematics.trajectories);
fprintf('Loaded. digitisation: 1x%d  |  trajectories: %d\n', numel(digitisation), nTrials);

% Trials to run
if isempty(opts.trialIndices)
    trialList = 1:nTrials;
else
    trialList = opts.trialIndices(:)';
end

% STL paths
femurStlPath = fullfile(projectRoot, 'model', 'cylinder', 'top_matlab1.stl');
tibiaStlPath = fullfile(projectRoot, 'model', 'cylinder', 'bottom_matlab1.stl');
assert(exist(femurStlPath, 'file') == 2, 'Cannot find top_matlab1.stl');
assert(exist(tibiaStlPath, 'file') == 2, 'Cannot find bottom_matlab1.stl');

% Read STLs (shared across all trials)
femur_mesh  = stlread(femurStlPath);
tibia_mesh  = stlread(tibiaStlPath);

V_femur_raw = femur_mesh.Points;   F_femur = femur_mesh.ConnectivityList;
V_tibia_raw = tibia_mesh.Points;   F_tibia = tibia_mesh.ConnectivityList;

fprintf('\nSTL info:\n');
fprintf('  Femur: %d vertices, X[%.1f %.1f] Y[%.1f %.1f] Z[%.1f %.1f]\n', ...
    size(V_femur_raw,1), min(V_femur_raw(:,1)), max(V_femur_raw(:,1)), ...
    min(V_femur_raw(:,2)), max(V_femur_raw(:,2)), min(V_femur_raw(:,3)), max(V_femur_raw(:,3)));
fprintf('  Tibia: %d vertices, X[%.1f %.1f] Y[%.1f %.1f] Z[%.1f %.1f]\n', ...
    size(V_tibia_raw,1), min(V_tibia_raw(:,1)), max(V_tibia_raw(:,1)), ...
    min(V_tibia_raw(:,2)), max(V_tibia_raw(:,2)), min(V_tibia_raw(:,3)), max(V_tibia_raw(:,3)));

% Apply STL errors (optional)
V_femur = V_femur_raw;
V_tibia = V_tibia_raw;

if ERR.active
    if ERR.stl_rotate_90_about_x
        Rx90 = [1 0 0; 0 0 -1; 0 1 0];
        V_femur = (Rx90 * V_femur')';
        fprintf('[ERROR INJECTION] Femur STL rotated 90deg about X\n');
    end
    if ERR.stl_flip_z
        V_femur(:,3) = -V_femur(:,3);
        fprintf('[ERROR INJECTION] Femur STL Z axis flipped\n');
    end
    if ERR.stl_offset_z_mm ~= 0
        V_femur(:,3) = V_femur(:,3) + ERR.stl_offset_z_mm;
        fprintf('[ERROR INJECTION] Femur STL shifted %.1f mm in Z\n', ERR.stl_offset_z_mm);
    end
end

% Metric names (as you had originally)
metrics = {'flexion',  'Flexion (deg)'; ...
           'varus',    'Varus (deg)'; ...
           'external', 'External rotation (deg)'; ...
           'lateral',  'Lateral (mm)'; ...
           'anterior', 'Anterior (mm)'; ...
           'superior', 'Superior (mm)'};

%% ====================== FIGURE 1: Digitisation =====================
if opts.doDigitisationFigure
    digitisation.visualise();
end

%% ================== FIGURE 2: Alignment Preview ====================
% Preview is a "frame 1 sanity check".
% - If animMode='tibiaFrame': show tibia fixed, femur transformed by tTf=inv(fTt)
% - If animMode='global'    : show BOTH transformed into global using gTti/gTfi
if opts.doAlignmentPreview
    pTrial = min(max(opts.previewTrial,1), nTrials);
    tr0 = kinematics.trajectories(pTrial);

    gTfi0 = tr0.Transform.gTfi;
    gTti0 = tr0.Transform.gTti;
    fTt0  = tr0.Transform.fTt;

    if ERR.active && ERR.fTt_swap_tibia_femur
        tmp = gTfi0; gTfi0 = gTti0; gTti0 = tmp;
        fprintf('[ERROR INJECTION] (preview) gTfi and gTti swapped\n');
    end

    figure('Name','Alignment preview (frame 1)','Color','w');
    hold on; grid on; axis equal; view(opts.viewAzEl(1), opts.viewAzEl(2));
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

    switch lower(opts.animMode)
        case 'tibiaframe'
            % tibia fixed in tibia frame; femur transformed to tibia frame
            tTf1 = inv(fTt0(:,:,1));
            patch('Faces', F_tibia, 'Vertices', V_tibia, ...
                'FaceColor', opts.tibiaColor, 'EdgeColor','none', 'FaceAlpha', 0.3);
            patch('Faces', F_femur, 'Vertices', applyRigid(tTf1, V_femur), ...
                'FaceColor', opts.femurColor, 'EdgeColor','none', 'FaceAlpha', 0.3);
            title(sprintf('Preview: tibia frame (trial %d, frame 1)', pTrial));

        case 'global'
            % both move into global
            Vt_g = applyRigid(gTti0(:,:,1), V_tibia);
            Vf_g = applyRigid(gTfi0(:,:,1), V_femur);

            patch('Faces', F_tibia, 'Vertices', Vt_g, ...
                'FaceColor', opts.tibiaColor, 'EdgeColor','none', 'FaceAlpha', 0.3);
            patch('Faces', F_femur, 'Vertices', Vf_g, ...
                'FaceColor', opts.femurColor, 'EdgeColor','none', 'FaceAlpha', 0.3);
            title(sprintf('Preview: global frame (trial %d, frame 1)', pTrial));

        otherwise
            error('opts.animMode must be ''tibiaFrame'' or ''global''.');
    end

    camlight; lighting gouraud;
    hold off;
end

%% ========================= LOOP OVER TRIALS =========================
for TRIAL_INDEX = trialList

    tr      = kinematics.trajectories(TRIAL_INDEX);
    gTfi    = tr.Transform.gTfi;
    gTti    = tr.Transform.gTti;
    fTt     = tr.Transform.fTt;

    nFrames = size(fTt,3);
    cond    = tr.LoadingCondition;

    fprintf('\n── Trial %d / %d : %s ──\n', TRIAL_INDEX, nTrials, cond);

    % ---- Apply transform errors (optional) ----
    if ERR.active
        if ERR.fTt_swap_tibia_femur
            tmp = gTfi; gTfi = gTti; gTti = tmp;
            fprintf('[ERROR INJECTION] gTfi and gTti swapped\n');
        end
        if ERR.fTt_add_flexion_deg ~= 0
            th = ERR.fTt_add_flexion_deg;
            Rx = [1 0 0 0; 0 cosd(th) -sind(th) 0; 0 sind(th) cosd(th) 0; 0 0 0 1];
            for k = 1:nFrames
                fTt(:,:,k) = Rx * fTt(:,:,k);
            end
            fprintf('[ERROR INJECTION] %.1f deg flexion added to fTt every frame\n', th);
        end
    end

    %% ---------- Body frame axes animation (global) ----------
    if opts.doBodyAxesAnimation
        femur_origins = squeeze(gTfi(1:3,4,:));
        tibia_origins = squeeze(gTti(1:3,4,:));

        SCALE = 30;  STEP = opts.stepAxes;
        figAxes = figure('Name', sprintf('Body frame axes — %s', cond), 'Color','w');

        for k = 1:STEP:nFrames
            clf(figAxes);
            hold on; grid on; axis equal; view(opts.viewAzEl(1), opts.viewAzEl(2));
            xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
            title(sprintf('%s — frame %d / %d', cond, k, nFrames));

            Rf = gTfi(1:3,1:3,k);  of = gTfi(1:3,4,k);
            Rt = gTti(1:3,1:3,k);  ot = gTti(1:3,4,k);

            for ax = 1:3
                quiver3(of(1),of(2),of(3), SCALE*Rf(1,ax),SCALE*Rf(2,ax),SCALE*Rf(3,ax), 0,'r','LineWidth',1.5);
                quiver3(ot(1),ot(2),ot(3), SCALE*Rt(1,ax),SCALE*Rt(2,ax),SCALE*Rt(3,ax), 0,'b','LineWidth',1.5);
            end
            text(of(1),of(2),of(3)+SCALE,'Femur','Color','r','FontWeight','bold');
            text(ot(1),ot(2),ot(3)+SCALE,'Tibia','Color','b','FontWeight','bold');

            plot3(femur_origins(1,1:k),femur_origins(2,1:k),femur_origins(3,1:k),'r:','LineWidth',0.5);
            plot3(tibia_origins(1,1:k),tibia_origins(2,1:k),tibia_origins(3,1:k),'b:','LineWidth',0.5);

            drawnow; pause(opts.pauseAxes);
            hold off;
        end
    end

    %% ---------- Kinematics (all six metrics) ----------
    if opts.doKinematicsPlots
        bones.femur.i      = squeeze(gTfi(1:3,1,:));
        bones.femur.j      = squeeze(gTfi(1:3,2,:));
        bones.femur.k      = squeeze(gTfi(1:3,3,:));
        bones.femur.origin = squeeze(gTfi(1:3,4,:));
        bones.tibia.i      = squeeze(gTti(1:3,1,:));
        bones.tibia.j      = squeeze(gTti(1:3,2,:));
        bones.tibia.k      = squeeze(gTti(1:3,3,:));
        bones.tibia.origin = squeeze(gTti(1:3,4,:));

        [tf, ~] = Knee.grood_and_suntay(bones.femur, bones.tibia, [], fTt, [], config.is_right_knee);

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

    %% ---------- STL animation ----------
    if opts.doSTLAnimation

        figStl = figure('Name', sprintf('STL — %s', cond), 'Color', 'w');
        hold on; grid on; axis equal; view(opts.viewAzEl(1), opts.viewAzEl(2));
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

        % Precompute frame transforms depending on mode
        V_femur_h = [V_femur' ; ones(1, size(V_femur,1))];
        V_tibia_h = [V_tibia' ; ones(1, size(V_tibia,1))];

        switch lower(opts.animMode)
            case 'tibiaframe'
                % Femur in tibia frame
                tTf = pageinv(fTt);                       % [4x4xN]
                Vf_all = pagemtimes(tTf, V_femur_h);      % [4 x Nv_f x N]

                % Tibia fixed (already tibia frame)
                Vt0 = V_tibia;

                % Axis ranges
                all_femur = reshape(Vf_all(1:3,:,:), 3, [])';
                all_verts = [all_femur; Vt0];
                x_range = [min(all_verts(:,1)) max(all_verts(:,1))];
                y_range = [min(all_verts(:,2)) max(all_verts(:,2))];
                z_range = [min(all_verts(:,3)) max(all_verts(:,3))];

                hTibia = patch('Faces',F_tibia,'Vertices',Vt0, ...
                    'FaceColor',opts.tibiaColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

                hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)', ...
                    'FaceColor',opts.femurColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

            case 'global'
                % Both in global
                Vf_all = pagemtimes(gTfi, V_femur_h);     % [4 x Nv_f x N]
                Vt_all = pagemtimes(gTti, V_tibia_h);     % [4 x Nv_t x N]

                all_femur = reshape(Vf_all(1:3,:,:), 3, [])';
                all_tibia = reshape(Vt_all(1:3,:,:), 3, [])';
                all_verts = [all_femur; all_tibia];

                x_range = [min(all_verts(:,1)) max(all_verts(:,1))];
                y_range = [min(all_verts(:,2)) max(all_verts(:,2))];
                z_range = [min(all_verts(:,3)) max(all_verts(:,3))];

                hTibia = patch('Faces',F_tibia,'Vertices',Vt_all(1:3,:,1)', ...
                    'FaceColor',opts.tibiaColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

                hFemur = patch('Faces',F_femur,'Vertices',Vf_all(1:3,:,1)', ...
                    'FaceColor',opts.femurColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);

            otherwise
                error('opts.animMode must be ''tibiaFrame'' or ''global''.');
        end

        xlim(x_range); ylim(y_range); zlim(z_range);
        camlight; lighting gouraud;
        legend([hFemur hTibia],{'Femur','Tibia'});

        % Optional video writer (OFF by default)
        if opts.saveVideo
            videoFile = fullfile(projectRoot, sprintf('cylinder_%s.mp4', cond));
            vw = VideoWriter(videoFile, 'MPEG-4');
            vw.FrameRate = opts.frameRate;
            open(vw);
        end

        for k = 1:nFrames
            set(hFemur, 'Vertices', Vf_all(1:3,:,k)');

            if strcmpi(opts.animMode,'global')
                set(hTibia, 'Vertices', Vt_all(1:3,:,k)');
            end

            title(sprintf('%s — frame %d / %d', cond, k, nFrames));
            drawnow;

            if opts.saveVideo
                writeVideo(vw, getframe(figStl));
            end
        end

        if opts.saveVideo
            close(vw);
            fprintf('Saved: %s\n', videoFile);
        else
            fprintf('(No video saved for %s)\n', cond);
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