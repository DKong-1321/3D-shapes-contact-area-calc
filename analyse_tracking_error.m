% analyse_tracking_error.m
% -------------------------------------------------------------------------
% Optical-tracking error analysis for the "one solid block, two STL halves"
% validation experiment.
%
% Setup:
%   * One solid block glued in place with two trackers attached (one to
%     the "top" half, one to the "bottom" half).
%   * The two halves are digitally separated as two STLs that touch along
%     the slicing plane.
%   * The camera was moved around the working volume while the block
%     stayed completely still.
%
% Because the block did not physically move, the TRUE contact area is
% constant (the area of the slicing plane). Any frame-to-frame change
% in the computed contact area is pure tracking error propagated into
% the contact-area metric.
%
% FULL-FRAME VERSION: processes every frameStep-th frame of each
% trajectory. With frameStep = 4 -> 88 of 352 and 109 of 435 frames
% (~197 contact computations total). Reduce frameStep to 1 for every
% frame, increase for faster runs.
%
% Pipeline:
%   1. Load .mat and STLs
%   2. STL position diagnostic (3 figures, console output)
%   3. Per-frame contact-area loop for each trajectory with live area
%      trace + MP4 export
%
% Outputs:
%   * Fig 1   Raw STLs (as authored on disk)
%   * Fig 2   STLs after frame-1 transform
%   * Fig 3   Side-by-side
%   * Fig 4+  One animation figure per trajectory + clean area-trace PNG
%   * MP4 per trajectory
%   * CSV of per-frame contact area
%   * Console summary: mean A, std A, RMS deviation
% -------------------------------------------------------------------------

clear; clc; close all;

%% =========================== USER SETTINGS ============================
primaryFile  = 'matlab 2.mat';
projectRoot  = fileparts(mfilename('fullpath'));

% STL paths (two halves of the same solid block)
femurStlPath = 'Femur-Top 2.stl';
tibiaStlPath = 'Tibia-bottom 2.stl';

% Contact function and tolerance
contactFcn    = @computeContactArea_STS_hybrid;
contactTol_mm = 0.50;

% Full-frame run
testRun     = false;         % full pass
TEST_FRAMES = 10;            % unused when testRun = false
frameStep   = 4;             % 1 = every frame, higher = sparser

% Visual / video settings
opts.viewAzEl     = [30 20];
opts.faceAlpha    = 0.85;
opts.topColor     = [0.85 0.70 0.55];
opts.botColor     = [0.55 0.70 0.85];
opts.contactColor = [0.90 0.20 0.20];
opts.saveVideo    = true;
opts.frameRate    = 8;

% STL position diagnostic
opts.runStlDiagnostic = true;

%% =============================== LOAD =================================
fprintf('\n=== analyse_tracking_error (full-frame) ===\n');
assert(exist(primaryFile,'file')==2, 'Cannot find %s', primaryFile);
S = load(primaryFile);
K = S.kinematics;
fprintf('Loaded %s\n', primaryFile);

try
    addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'lib')));
    addpath(genpath(fullfile(projectRoot, 'OpticalTracking-dev-3', 'external')));
catch
end

fcnName = func2str(contactFcn);
assert(exist(fcnName,'file')==2, ...
    'Contact function "%s" not on MATLAB path. addpath the folder containing it.', fcnName);

%% =================== STL LOADING ======================================
[F_top, V_top] = tryLoadStl(femurStlPath);
[F_bot, V_bot] = tryLoadStl(tibiaStlPath);
assert(~isempty(F_top), 'Cannot load top STL: %s', femurStlPath);
assert(~isempty(F_bot), 'Cannot load bottom STL: %s', tibiaStlPath);
fprintf('Loaded STLs:\n');
fprintf('  Top    (%s): %d verts / %d faces\n', femurStlPath, size(V_top,1), size(F_top,1));
fprintf('  Bottom (%s): %d verts / %d faces\n', tibiaStlPath, size(V_bot,1), size(F_bot,1));

%% ==================== STL POSITION DIAGNOSTIC =========================
if opts.runStlDiagnostic
    runStlPositionDiagnostic(K, F_top, V_top, F_bot, V_bot, opts);
end

%% =================== PROCESS BOTH TRAJECTORIES ========================
trajList = K.trajectories;

% Estimate runtime up front so you know what you're committing to
totalFrames = 0;
for t = 1:numel(trajList)
    nF = size(trajList(t).Transform.gTfi, 3);
    totalFrames = totalFrames + numel(1:frameStep:nF);
end
fprintf('\nFull-frame run plan:\n');
fprintf('  frameStep = %d  ->  ~%d total contact computations across all trajectories\n', ...
    frameStep, totalFrames);
fprintf('  Rough estimate: with 2 s per contact call, expect ~%.1f min total.\n', ...
    2*totalFrames/60);
fprintf('  (Press Ctrl+C to abort; partial results are not saved until each trajectory finishes.)\n');

for t = 1:numel(trajList)
    name    = char(trajList(t).LoadingCondition);
    gTfi    = trajList(t).Transform.gTfi;
    gTti    = trajList(t).Transform.gTti;
    nFrames = size(gTfi,3);

    frames = 1:frameStep:nFrames;
    fprintf('\n=== Trajectory %d: %s ===\n   Every %dth frame: %d of %d\n', ...
        t, name, frameStep, numel(frames), nFrames);

    runOneTrajectory(t, name, gTfi, gTti, frames, ...
        F_top, V_top, F_bot, V_bot, ...
        contactFcn, contactTol_mm, opts, projectRoot);
end

fprintf('\nDone.\n');

%% =======================================================================
%% =========================== HELPERS ===================================
%% =======================================================================
function [F, V] = tryLoadStl(p)
F = []; V = [];
if isempty(p) || exist(p,'file') ~= 2, return; end
try
    m = stlread(p);
    F = m.ConnectivityList;
    V = m.Points;
catch
end
end

% -------------------------------------------------------------------------
function Vout = applyRigid(T, V)
Vh = [V ones(size(V,1),1)];
Vt = (T * Vh.').';
Vout = Vt(:,1:3);
end

% -------------------------------------------------------------------------
function runStlPositionDiagnostic(K, F_top, V_top, F_bot, V_bot, opts)
tr1   = K.trajectories(1);
gTfi1 = tr1.Transform.gTfi(:,:,1);
gTti1 = tr1.Transform.gTti(:,:,1);

fprintf('\nSTL POSITION DIAGNOSTIC\n');
fprintf('-----------------------\n');
fprintf('Raw STLs (as authored on disk):\n');
printStlInfo('  Top    half', V_top);
printStlInfo('  Bottom half', V_bot);
gapRaw = mean(V_top,1) - mean(V_bot,1);
fprintf('  Raw centroid gap (top - bottom) = [%7.2f %7.2f %7.2f]  |gap| = %.2f mm\n', ...
    gapRaw, norm(gapRaw));

V_top_g = applyRigid(gTfi1, V_top);
V_bot_g = applyRigid(gTti1, V_bot);

fprintf('\nTransforms at frame 1 of trajectory 1:\n');
fprintf('  gTfi(:,:,1)  (top  half):\n');  disp(gTfi1);
fprintf('  gTti(:,:,1)  (bottom half):\n'); disp(gTti1);

fprintf('After frame-1 transform (global coords):\n');
printStlInfo('  Top    half', V_top_g);
printStlInfo('  Bottom half', V_bot_g);
gapG = mean(V_top_g,1) - mean(V_bot_g,1);
fprintf('  Transformed centroid gap = [%7.2f %7.2f %7.2f]  |gap| = %.2f mm\n', ...
    gapG, norm(gapG));

sample = randperm(size(V_top_g,1), min(500,size(V_top_g,1)));
kdb = KDTreeSearcher(V_bot_g);
[~, d_t2b] = knnsearch(kdb, V_top_g(sample,:));
fprintf('\n  Closest-point distance top->bottom: min=%.3f  median=%.3f  max=%.3f mm\n', ...
    min(d_t2b), median(d_t2b), max(d_t2b));

% Fig 1: raw STLs
figure('Color','w','Name','Fig 1 - Raw STLs','Position',[60 60 700 600]);
hold on; grid on; axis equal; view(opts.viewAzEl);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
patch('Faces',F_top,'Vertices',V_top,'FaceColor',opts.topColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
patch('Faces',F_bot,'Vertices',V_bot,'FaceColor',opts.botColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
plotBBox(V_top, opts.topColor);
plotBBox(V_bot, opts.botColor);
camlight; lighting gouraud;
title({'Fig 1: Raw STL files as authored on disk', ...
    sprintf('raw centroid gap = %.2f mm', norm(gapRaw))});

% Fig 2: frame-1 transformed
figure('Color','w','Name','Fig 2 - Frame 1 transformed','Position',[60 60 700 600]);
hold on; grid on; axis equal; view(opts.viewAzEl);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
patch('Faces',F_top,'Vertices',V_top_g,'FaceColor',opts.topColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
patch('Faces',F_bot,'Vertices',V_bot_g,'FaceColor',opts.botColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
plotBBox(V_top_g, opts.topColor);
plotBBox(V_bot_g, opts.botColor);
camlight; lighting gouraud;
title({'Fig 2: STLs after frame-1 transform', ...
    sprintf('transformed gap = %.2f mm, closest-pt min = %.3f mm', ...
    norm(gapG), min(d_t2b))});

% Fig 3: side-by-side
figure('Color','w','Name','Fig 3 - Side-by-side','Position',[60 60 1300 600]);
subplot(1,2,1); hold on; grid on; axis equal; view(opts.viewAzEl);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
patch('Faces',F_top,'Vertices',V_top,'FaceColor',opts.topColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
patch('Faces',F_bot,'Vertices',V_bot,'FaceColor',opts.botColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
camlight; lighting gouraud;
title(sprintf('Raw STLs  (gap %.2f mm)', norm(gapRaw)));

subplot(1,2,2); hold on; grid on; axis equal; view(opts.viewAzEl);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
patch('Faces',F_top,'Vertices',V_top_g,'FaceColor',opts.topColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
patch('Faces',F_bot,'Vertices',V_bot_g,'FaceColor',opts.botColor, ...
    'EdgeColor','none','FaceAlpha',0.55);
camlight; lighting gouraud;
title(sprintf('Frame-1 transformed  (gap %.2f mm)', norm(gapG)));

sgtitle('Left = raw STL files, Right = after frame-1 transform');
end

% -------------------------------------------------------------------------
function printStlInfo(label, V)
bbox   = [min(V); max(V)];
extent = bbox(2,:) - bbox(1,:);
c = mean(V,1);
fprintf('%s : %5d verts  bbox X[%7.2f %7.2f] Y[%7.2f %7.2f] Z[%7.2f %7.2f]\n', ...
    label, size(V,1), bbox(1,1),bbox(2,1), bbox(1,2),bbox(2,2), bbox(1,3),bbox(2,3));
fprintf('%s   extent [%6.2f %6.2f %6.2f] mm  centroid [%7.2f %7.2f %7.2f]\n', ...
    repmat(' ',1,numel(label)), extent(1),extent(2),extent(3), c(1),c(2),c(3));
end

% -------------------------------------------------------------------------
function plotBBox(V, col)
mn = min(V); mx = max(V);
corners = [mn(1) mn(2) mn(3); mx(1) mn(2) mn(3); mx(1) mx(2) mn(3); mn(1) mx(2) mn(3); ...
           mn(1) mn(2) mx(3); mx(1) mn(2) mx(3); mx(1) mx(2) mx(3); mn(1) mx(2) mx(3)];
edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
for i = 1:size(edges,1)
    plot3(corners(edges(i,:),1), corners(edges(i,:),2), corners(edges(i,:),3), ...
        '-', 'Color', col, 'LineWidth', 0.6);
end
end

% -------------------------------------------------------------------------
function runOneTrajectory(tIdx, condName, gTfi, gTti, frames, ...
        F_top, V_top, F_bot, V_bot, ...
        contactFcn, tol_mm, opts, projectRoot)

nF = numel(frames);

V_top_h = [V_top, ones(size(V_top,1),1)].';
V_bot_h = [V_bot, ones(size(V_bot,1),1)].';
Vtop_all = pagemtimes(gTfi(:,:,frames), V_top_h);
Vbot_all = pagemtimes(gTti(:,:,frames), V_bot_h);

allV = [reshape(Vtop_all(1:3,:,:),3,[]).'; reshape(Vbot_all(1:3,:,:),3,[]).'];
pad  = 5;
xlims = [min(allV(:,1))-pad max(allV(:,1))+pad];
ylims = [min(allV(:,2))-pad max(allV(:,2))+pad];
zlims = [min(allV(:,3))-pad max(allV(:,3))+pad];

fig = figure('Color','w','Position',[60 60 1200 820], ...
    'Name',sprintf('Tracking-error contact area - %s', condName));
ax3 = subplot(3,1,[1 2],'Parent',fig);
hold(ax3,'on'); grid(ax3,'on'); axis(ax3,'equal');
view(ax3, opts.viewAzEl(1), opts.viewAzEl(2));
xlabel(ax3,'X (mm)'); ylabel(ax3,'Y (mm)'); zlabel(ax3,'Z (mm)');
xlim(ax3,xlims); ylim(ax3,ylims); zlim(ax3,zlims);

hTop = patch(ax3,'Faces',F_top,'Vertices',Vtop_all(1:3,:,1).', ...
    'FaceColor',opts.topColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);
hBot = patch(ax3,'Faces',F_bot,'Vertices',Vbot_all(1:3,:,1).', ...
    'FaceColor',opts.botColor,'EdgeColor','none','FaceAlpha',opts.faceAlpha);
hCon = patch(ax3,'Faces', zeros(0,3),'Vertices',Vbot_all(1:3,:,1).', ...
    'FaceColor',opts.contactColor,'EdgeColor','none','FaceAlpha',1.0);

camlight(ax3); lighting(ax3,'gouraud');
legend(ax3, [hTop hBot hCon], {'Top half','Bottom half','Contact patch'}, ...
    'Location','northeast');

axA = subplot(3,1,3,'Parent',fig);
hold(axA,'on'); grid(axA,'on');
xlabel(axA,'Frame'); ylabel(axA,'Contact area (mm^2)');
hLine = plot(axA, nan, nan, '-o', 'LineWidth',1.4, ...
    'Color',[0.20 0.40 0.85], 'MarkerFaceColor',[0.20 0.40 0.85], 'MarkerSize',5);
hMark = plot(axA, nan, nan, 'ro', 'MarkerSize',9, 'LineWidth',1.5);
xlim(axA, [frames(1)-1 frames(end)+1]);

hTitle = sgtitle(fig, sprintf('%s   |   tol = %.3f mm   |   frame %d / %d', ...
    condName, tol_mm, frames(1), frames(end)), 'Interpreter','none');

vw = [];
if opts.saveVideo
    safeName  = matlab.lang.makeValidName(condName);
    videoFile = fullfile(projectRoot, sprintf('tracking_contact_%s.mp4', safeName));
    try, vw = VideoWriter(videoFile, 'MPEG-4');
    catch
        videoFile = fullfile(projectRoot, sprintf('tracking_contact_%s.avi', safeName));
        vw = VideoWriter(videoFile);
    end
    vw.FrameRate = opts.frameRate;
    open(vw);
end

A_per_frame = nan(nF,1);
t0 = tic;

for k = 1:nF
    V_top_k = Vtop_all(1:3,:,k).';
    V_bot_k = Vbot_all(1:3,:,k).';

    master = buildBodyStruct(F_top, V_top_k);
    slave  = buildBodyStruct(F_bot, V_bot_k);

    res = [];
    try
        res = contactFcn(master, slave, tol_mm);
    catch ME1
        try
            res = contactFcn(master, slave, struct('tol', tol_mm));
        catch ME2
            error(['Contact function call failed both ways:\n' ...
                   '  (master,slave,tol)  -> %s\n' ...
                   '  (master,slave,opts) -> %s\n'], ME1.message, ME2.message);
        end
    end

    A_k = pickField(res, {'contactArea','area','A_total','A','contact_area'});
    if isnan(A_k)
        error('No area field on result. Fields: %s', strjoin(fieldnames(res), ', '));
    end
    A_per_frame(k) = A_k;

    mask = pickField(res, {'contactMask','mask','contact_mask','inContact','slaveMask'});

    set(hTop, 'Vertices', V_top_k);
    set(hBot, 'Vertices', V_bot_k);
    if isnumeric(mask) && ~all(isnan(mask(:))) && numel(mask) == size(F_bot,1)
        set(hCon, 'Faces', F_bot(logical(mask),:), 'Vertices', V_bot_k);
    elseif islogical(mask) && numel(mask) == size(F_bot,1)
        set(hCon, 'Faces', F_bot(mask,:), 'Vertices', V_bot_k);
    else
        set(hCon, 'Faces', zeros(0,3));
    end

    set(hLine, 'XData', frames(1:k), 'YData', A_per_frame(1:k));
    set(hMark, 'XData', frames(k),   'YData', A_per_frame(k));
    set(hTitle,'String', sprintf('%s   |   tol=%.3f mm   |   frame %d/%d   |   A = %.2f mm^2', ...
        condName, tol_mm, frames(k), frames(end), A_per_frame(k)));
    drawnow;

    if ~isempty(vw), writeVideo(vw, getframe(fig)); end

    % Console: print every frame in test runs, less often for full runs
    if nF <= 20 || mod(k,10)==0 || k == nF
        elapsed = toc(t0);
        eta = elapsed * (nF - k) / max(k,1);
        fprintf('  frame %4d / %4d   (recording frame %4d)   A = %8.2f mm^2   elapsed=%5.1fs  ETA=%5.1fs\n', ...
            k, nF, frames(k), A_per_frame(k), elapsed, eta);
    end
end

if ~isempty(vw)
    close(vw);
    fprintf('  Saved video: tracking_contact_%s.mp4\n', condName);
end

A_mean = mean(A_per_frame,'omitnan');
A_std  = std(A_per_frame,0,'omitnan');
A_rms  = sqrt(mean((A_per_frame - A_mean).^2,'omitnan'));
A_p95  = prctile(abs(A_per_frame - A_mean), 95);
A_min  = min(A_per_frame,[],'omitnan');
A_max  = max(A_per_frame,[],'omitnan');

fprintf('\nTrajectory %d (%s) - contact-area summary (tol = %.3f mm):\n', tIdx, condName, tol_mm);
fprintf('  mean A          = %.3f mm^2\n', A_mean);
fprintf('  std A           = %.3f mm^2   (tracking error in mm^2)\n', A_std);
fprintf('  RMS deviation   = %.3f mm^2\n', A_rms);
fprintf('  95th pct |dev|  = %.3f mm^2\n', A_p95);
fprintf('  range [min max] = [%.3f  %.3f] mm^2   (max - min = %.3f)\n', A_min, A_max, A_max - A_min);
fprintf('  rel std         = %.2f %%   (std / mean)\n', 100 * A_std / A_mean);

safeName = matlab.lang.makeValidName(condName);
csvFile = fullfile(projectRoot, sprintf('contact_area_per_frame_%s.csv', safeName));
Tcsv = table(frames(:), A_per_frame, 'VariableNames', {'frame','contactArea_mm2'});
writetable(Tcsv, csvFile);
fprintf('  Saved CSV   : %s\n', csvFile);

figA = figure('Color','w','Name',sprintf('Contact area trace - %s', condName));
plot(frames, A_per_frame, '-o', 'LineWidth',1.4); grid on; hold on;
yline(A_mean,'--k','LineWidth',1.2,'DisplayName',sprintf('mean = %.2f mm^2', A_mean));
xlabel('Frame'); ylabel('Contact area (mm^2)');
title(sprintf('%s - contact area vs frame  (mean %.2f, std %.3f mm^2)', ...
    condName, A_mean, A_std), 'Interpreter','none');
legend('Location','best');
pngFile = fullfile(projectRoot, sprintf('contact_area_trace_%s.png', safeName));
exportgraphics(figA, pngFile, 'Resolution', 200);
fprintf('  Saved PNG   : %s\n', pngFile);
end

% -------------------------------------------------------------------------
function v = pickField(s, candidates)
v = nan;
for i = 1:numel(candidates)
    if isfield(s, candidates{i})
        v = s.(candidates{i}); return;
    end
end
end