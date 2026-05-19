% analyse_errors.m
% -------------------------------------------------------------------------
% Final error-budget analysis for the optical-tracking contact-area
% pipeline. Pulls every source of error from a single .mat file :
%
%   1) WITHIN-TOUCH SENSOR NOISE   — spread of the 18 Tx/Ty/Tz samples
%      of each landmark held by the probe fzor ~1 s. This is the noise
%      of the optical system PLUS the probe operator's micro-tremor.
%
%   2) BETWEEN-TOUCH DIGITISATION ERROR — how much the SAME landmark
%      moves between separately-recorded touch sessions (the tl1, tl2,
%      tl3, tl4, tl5 mentioned). Needs the other 4 .mat files.
%      Leave digitisationFiles = {} to skip this.
%
%   3) OPTICAL-TRACKING ERROR — body clamped, camera moved around its
%      working limits. Lives inside
%         kinematics.trajectories(k).Transform.gTfi  (4 x 4 x Nframes)
%         kinematics.trajectories(k).Transform.gTti
%      Any frame-to-frame variation in those poses, with the body still,
%      is the camera's own tracking error.
%
%   4) ERROR vs TOLERANCE + analytical solution
%      Loads toleranceSweep_*.csv (which already has A_true_mm2 per pair)
%      and overlays the analytical area on each numerical curve, with an
%      error band derived from the combined sigma above.
%
% Output:
%   Fig 1  per-landmark within-touch RMS noise
%   Fig 2  per-landmark between-session digitisation RMS (if N_sessions>=2)
%   Fig 3  per-trajectory tracking error: bone translation & rotation drift
%          versus frame index (and overall RMS)
%   Fig 4  summary boxplot of all three sources side by side
%   Fig 5  contact-area vs tolerance with analytical A_true + sigma band
%   error_summary_<timestamp>.csv
% -------------------------------------------------------------------------

clear; clc; close all;

%% =========================== USER SETTINGS ============================
primaryFile        = 'matlab 2.mat';                     
digitisationFiles  = {};                                 % add the other 4 .mat paths here later
csvFile            = 'toleranceSweep_20260126_142129.csv'; % '' to skip Fig 5
excludeLandmarks   = {'tracker'};                        % skip rigid-body markers
bonesToUse         = {'tibia','femur'};                  % drop patella unless needed

%% =============================== LOAD =================================
fprintf('\n=== analyse_errors ===\n');
assert(exist(primaryFile,'file')==2, 'Cannot find %s', primaryFile);
S  = load(primaryFile);
D  = S.digitisation;
K  = S.kinematics;
fprintf('Loaded %s\n', primaryFile);

%% ============ DISCOVER LANDMARKS, DROP TRACKERS/PATELLA ===============
landmarkKeys = listLandmarks(D.bone, bonesToUse, excludeLandmarks);
fprintf('Using %d landmarks: %s\n', numel(landmarkKeys), strjoin(landmarkKeys,', '));

%% =================== 1) WITHIN-TOUCH SENSOR NOISE =====================
fprintf('\n--- (1) Within-touch sensor noise (18-sample hold per landmark) ---\n');
withinStats = withinTouchStats(D.bone, landmarkKeys);
printStats(withinStats);

%% =================== 2) BETWEEN-TOUCH DIGITISATION ====================
allFiles = [{primaryFile} digitisationFiles(:)'];
nSessions = numel(allFiles);
sessionMeans = nan(numel(landmarkKeys), 3, nSessions);
for s = 1:nSessions
    fprintf('Reading session %d/%d: %s\n', s, nSessions, allFiles{s});
    Sj = load(allFiles{s});
    Bj = Sj.digitisation.bone;
    for k = 1:numel(landmarkKeys)
        sessionMeans(k,:,s) = touchMean(Bj, landmarkKeys{k});
    end
end

if nSessions >= 2
    fprintf('\n--- (2) Between-session digitisation error (across %d sessions) ---\n', nSessions);
    digStats = perLandmarkStatsFromStack(sessionMeans, landmarkKeys);
    printStats(digStats);
else
    fprintf('\n--- (2) Between-session digitisation error: skipped (need >=2 .mat files) ---\n');
    digStats = emptyStats(landmarkKeys);
end

%% =================== 3) OPTICAL-TRACKING ERROR ========================
fprintf('\n--- (3) Optical-tracking error (model still, camera moved) ---\n');
trajList = K.trajectories;
trackResults = struct([]);
for t = 1:numel(trajList)
    name = string(trajList(t).LoadingCondition);
    fprintf('\nTrajectory %d: %s\n', t, name);
    Tr = trajList(t).Transform;
    [stat_f, posF, rotF] = trackingErrorFromHomogeneous(Tr.gTfi, 'femur');
    [stat_t, posT, rotT] = trackingErrorFromHomogeneous(Tr.gTti, 'tibia');
    trackResults(t).name      = name;
    trackResults(t).femur     = stat_f;
    trackResults(t).tibia     = stat_t;
    trackResults(t).femur_pos = posF;   trackResults(t).femur_rot = rotF;
    trackResults(t).tibia_pos = posT;   trackResults(t).tibia_rot = rotT;
end

% Combined tracking sigma: RMS of all per-frame position residuals across
% all bones and trajectories (rotations reported separately for context).
allPosRes = [];
for t = 1:numel(trackResults)
    allPosRes = [allPosRes; trackResults(t).femur.posResMag(:); ...
                            trackResults(t).tibia.posResMag(:)]; %#ok<AGROW>
end
sigma_track = sqrt(mean(allPosRes.^2,'omitnan'));

%% =================== COMBINED BUDGET ==================================
sigma_within = mean(withinStats.rms_mm,'omitnan');
sigma_dig    = mean(digStats.rms_mm,   'omitnan');
if isnan(sigma_dig), sigma_dig = sigma_within; end   % fall back to within-touch
sigma_total  = sqrt(sigma_dig^2 + sigma_track^2);

fprintf('\n=== Combined error budget ===\n');
fprintf('  sigma_within-touch sensor noise : %.4f mm  (n_landmarks)\n', sigma_within);
fprintf('  sigma_between-session digitis.  : %.4f mm  (%s)\n', sigma_dig, ...
    ternary(nSessions>=2, sprintf('%d sessions', nSessions), 'fallback = within-touch'));
fprintf('  sigma_optical-tracking          : %.4f mm  (n=%d pose frames)\n', sigma_track, numel(allPosRes));
fprintf('  sigma_combined = sqrt(d^2 + t^2): %.4f mm\n', sigma_total);

%% =================== FIGURES ==========================================
plotPerLandmarkBars(withinStats, 'Within-touch sensor noise per landmark (18-sample hold)', 1, [0.20 0.40 0.85]);
if nSessions >= 2
    plotPerLandmarkBars(digStats, sprintf('Between-session digitisation error (%d sessions)', nSessions), 2, [0.30 0.65 0.30]);
end
plotTrackingError(trackResults);
plotSummaryBox(withinStats, digStats, trackResults);

if ~isempty(csvFile) && exist(csvFile,'file')==2
    plotToleranceSweep(csvFile, sigma_total);
else
    fprintf('\nNo tolerance CSV found at "%s" - skipping Fig 5.\n', csvFile);
end

%% =================== WRITE SUMMARY CSV ================================
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outCsv = sprintf('error_summary_%s.csv', ts);
writeSummary(outCsv, withinStats, digStats, trackResults, sigma_within, sigma_dig, sigma_track, sigma_total);
fprintf('\nSaved: %s\n', outCsv);

%% =======================================================================
%% =========================== HELPERS ===================================
%% =======================================================================
function keys = listLandmarks(bone, bonesToUse, excludeLandmarks)
keys = {};
bones = fieldnames(bone);
if ~isempty(bonesToUse), bones = intersect(bones, bonesToUse, 'stable'); end
for b = 1:numel(bones)
    sub = bone.(bones{b});
    if ~isstruct(sub) && ~isobject(sub), continue; end
    lms = fieldnames(sub);
    for l = 1:numel(lms)
        if any(strcmpi(lms{l}, excludeLandmarks)), continue; end
        keys{end+1} = sprintf('%s.%s', bones{b}, lms{l}); %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function tr = getTracker(bone, key)
parts = strsplit(key,'.');
obj = bone;
for i = 1:numel(parts)
    obj = obj.(parts{i});
end
tr = obj.unwrap();   % Tracker
end

% -------------------------------------------------------------------------
function mu = touchMean(bone, key)
mu = [nan nan nan];
try
    tr = getTracker(bone, key);
    P = [tr.Tx tr.Ty tr.Tz];
    P = P(all(~isnan(P),2),:);
    if ~isempty(P), mu = mean(P,1); end
catch
end
end

% -------------------------------------------------------------------------
function st = withinTouchStats(bone, keys)
N = numel(keys);
st = emptyStats(keys);
for n = 1:N
    try
        tr = getTracker(bone, keys{n});
        P = [tr.Tx tr.Ty tr.Tz];
        P = P(all(~isnan(P),2),:);
        if size(P,1) < 2, continue; end
        mu = mean(P,1); R = P - mu; d = sqrt(sum(R.^2,2));
        st.mean_xyz(n,:)   = mu;
        st.rms_mm(n)       = sqrt(mean(d.^2));
        st.p95_mm(n)       = prctile(d,95);
        st.maxRes_mm(n)    = max(d);
        st.std_xyz_mm(n,:) = std(P,0,1);
        st.nRepeats(n)     = size(P,1);
    catch
    end
end
end

% -------------------------------------------------------------------------
function st = perLandmarkStatsFromStack(stack, keys)
% stack: Nlm x 3 x Nsessions  (one mean per landmark per session)
N = numel(keys);
st = emptyStats(keys);
for n = 1:N
    P = squeeze(stack(n,:,:)).';       % Ns x 3
    P = P(all(~isnan(P),2),:);
    if size(P,1) < 2, continue; end
    mu = mean(P,1); R = P - mu; d = sqrt(sum(R.^2,2));
    st.mean_xyz(n,:)   = mu;
    st.rms_mm(n)       = sqrt(mean(d.^2));
    st.p95_mm(n)       = prctile(d,95);
    st.maxRes_mm(n)    = max(d);
    st.std_xyz_mm(n,:) = std(P,0,1);
    st.nRepeats(n)     = size(P,1);
end
end

% -------------------------------------------------------------------------
function [st, posResMag, rotResDeg] = trackingErrorFromHomogeneous(gT, label)
% gT: 4x4xN homogeneous transforms of one bone tracker, body still.
% Returns:
%   posResMag: Nx1 magnitudes of (position - mean position)
%   rotResDeg: Nx1 geodesic rotation angles from mean orientation, deg
N = size(gT,3);
pos = squeeze(gT(1:3,4,:)).';          % N x 3
muPos = mean(pos,1,'omitnan');
posRes = pos - muPos;
posResMag = sqrt(sum(posRes.^2,2));

% Mean orientation: simple averaging then re-orthogonalise via SVD
Rstack = gT(1:3,1:3,:);
Rmean = mean(Rstack,3,'omitnan');
[U,~,V] = svd(Rmean);
Rbar = U*V.';
if det(Rbar) < 0
    V(:,end) = -V(:,end); Rbar = U*V.';
end

rotResDeg = nan(N,1);
for i = 1:N
    Ri = Rstack(:,:,i);
    if any(isnan(Ri(:))), continue; end
    Rd = Rbar.' * Ri;
    cosang = (trace(Rd) - 1)/2;
    cosang = max(-1, min(1, cosang));
    rotResDeg(i) = rad2deg(acos(cosang));
end

st.label       = label;
st.posRms_mm   = sqrt(mean(posResMag.^2,'omitnan'));
st.posP95_mm   = prctile(posResMag,95);
st.posMax_mm   = max(posResMag,[],'omitnan');
st.rotRms_deg  = sqrt(mean(rotResDeg.^2,'omitnan'));
st.rotP95_deg  = prctile(rotResDeg,95);
st.rotMax_deg  = max(rotResDeg,[],'omitnan');
st.nFrames     = sum(~isnan(posResMag));
st.posResMag   = posResMag;
st.rotResDeg   = rotResDeg;
end

% -------------------------------------------------------------------------
function st = emptyStats(keys)
N = numel(keys);
st = struct('name',{keys(:)}, ...
            'mean_xyz',    nan(N,3), ...
            'rms_mm',      nan(N,1), ...
            'p95_mm',      nan(N,1), ...
            'maxRes_mm',   nan(N,1), ...
            'std_xyz_mm',  nan(N,3), ...
            'nRepeats',    nan(N,1));
end

% -------------------------------------------------------------------------
function printStats(st)
fprintf('%-22s %6s %10s %10s %10s\n','landmark','n','RMS','p95','max');
for n = 1:numel(st.name)
    if isnan(st.rms_mm(n)), continue; end
    fprintf('%-22s %6d %10.4f %10.4f %10.4f\n', ...
        st.name{n}, st.nRepeats(n), st.rms_mm(n), st.p95_mm(n), st.maxRes_mm(n));
end
fprintf('%-22s        %10.4f %10.4f %10.4f   (mean across landmarks)\n', ...
    'OVERALL', mean(st.rms_mm,'omitnan'), mean(st.p95_mm,'omitnan'), mean(st.maxRes_mm,'omitnan'));
end

% -------------------------------------------------------------------------
function plotPerLandmarkBars(st, ttl, figNum, color)
if all(isnan(st.rms_mm)), return; end
figure(figNum); clf; set(gcf,'Color','w');
bar(st.rms_mm, 'FaceColor', color, 'EdgeColor', 'none'); hold on;
for n = 1:numel(st.name)
    if isnan(st.p95_mm(n)), continue; end
    plot([n n], [st.rms_mm(n) st.p95_mm(n)], 'k-', 'LineWidth', 1);
    plot(n, st.p95_mm(n), 'k_', 'MarkerSize', 10, 'LineWidth', 1);
end
grid on;
set(gca,'XTick',1:numel(st.name),'XTickLabel',st.name, ...
    'XTickLabelRotation',45,'TickLabelInterpreter','none');
ylabel('Residual (mm)');
title({ttl,'(bar = RMS, cap = 95th percentile)'},'Interpreter','none');
end

% -------------------------------------------------------------------------
function plotTrackingError(R)
figure(3); clf; set(gcf,'Color','w');
nT = numel(R);
for t = 1:nT
    subplot(2, nT, t); hold on; grid on;
    plot(R(t).femur_pos,'-','LineWidth',1.2,'DisplayName','femur');
    plot(R(t).tibia_pos,'-','LineWidth',1.2,'DisplayName','tibia');
    xlabel('Frame'); ylabel('Position residual |\Delta r| (mm)');
    title({char(R(t).name), ...
        sprintf('femur RMS=%.3f mm,  tibia RMS=%.3f mm', R(t).femur.posRms_mm, R(t).tibia.posRms_mm)}, ...
        'Interpreter','none');
    legend('Location','best');

    subplot(2, nT, nT+t); hold on; grid on;
    plot(R(t).femur_rot,'-','LineWidth',1.2,'DisplayName','femur');
    plot(R(t).tibia_rot,'-','LineWidth',1.2,'DisplayName','tibia');
    xlabel('Frame'); ylabel('Rotation residual (deg)');
    title(sprintf('femur RMS=%.3f deg,  tibia RMS=%.3f deg', R(t).femur.rotRms_deg, R(t).tibia.rotRms_deg));
end
sgtitle('Optical-tracking error (body still, camera moved)');
end

% -------------------------------------------------------------------------
function plotSummaryBox(withinStats, digStats, trackResults)
data = []; group = {};
v = withinStats.rms_mm(~isnan(withinStats.rms_mm));
if ~isempty(v), data = [data; v]; group = [group; repmat({'Within-touch'},numel(v),1)]; end
v = digStats.rms_mm(~isnan(digStats.rms_mm));
if ~isempty(v), data = [data; v]; group = [group; repmat({'Between-session'},numel(v),1)]; end
allPos = [];
for t = 1:numel(trackResults)
    allPos = [allPos; trackResults(t).femur.posResMag(:); trackResults(t).tibia.posResMag(:)]; %#ok<AGROW>
end
allPos = allPos(~isnan(allPos));
if ~isempty(allPos), data = [data; allPos]; group = [group; repmat({'Tracking (per-frame)'},numel(allPos),1)]; end
if isempty(data), return; end
figure(4); clf; set(gcf,'Color','w');
boxplot(data, group); grid on;
set(gca,'YScale','log');
ylabel('Residual magnitude (mm)  [log scale]');
title('Error budget summary');
end

% -------------------------------------------------------------------------
function plotToleranceSweep(csvFile, sigma_total)
T = readtable(csvFile);
pairs = unique(string(T.pair),'stable');
figure(5); clf; set(gcf,'Color','w');
nP = numel(pairs);
for i = 1:nP
    subplot(nP,1,i); hold on; grid on;
    rows = string(T.pair) == pairs(i);
    tol  = T.tol_mm(rows); A = T.A_mm2(rows); A_tr = T.A_true_mm2(rows);
    [tol, idx] = sort(tol); A = A(idx); A_tr = A_tr(idx);

    plot(tol, A, '-o', 'LineWidth',1.5, 'Color',[0.20 0.40 0.85], ...
        'MarkerFaceColor',[0.20 0.40 0.85], 'DisplayName','Numerical A(tol)');
    yline(A_tr(1),'--k','LineWidth',1.2,'DisplayName',...
        sprintf('Analytical A_{true} = %.3f mm^2', A_tr(1)));
    dA = 2 * sqrt(pi * A_tr(1)) * sigma_total;
    yline(A_tr(1)+dA,':r','LineWidth',1.0,'DisplayName',...
        sprintf('\\pm 2\\surd(\\pi A)\\sigma  (\\sigma=%.3f mm)', sigma_total));
    yline(A_tr(1)-dA,':r','LineWidth',1.0,'HandleVisibility','off');

    xlabel('Contact tolerance (mm)'); ylabel('Contact area (mm^2)');
    title(sprintf('%s', pairs(i)), 'Interpreter','none');
    if i == 1, legend('Location','best'); end
end
sgtitle('Contact area vs tolerance - numerical, analytical, and propagated error band');
end

% -------------------------------------------------------------------------
function writeSummary(outCsv, withinStats, digStats, trackResults, sw, sd, sst, scomb)
rows = {};
    function addRows(label, stats)
        for n = 1:numel(stats.name)
            if isnan(stats.rms_mm(n)), continue; end
            rows(end+1,:) = { label, stats.name{n}, stats.nRepeats(n), ...
                stats.rms_mm(n), stats.p95_mm(n), stats.maxRes_mm(n) }; %#ok<AGROW>
        end
    end
addRows('within_touch', withinStats);
addRows('between_session', digStats);
for t = 1:numel(trackResults)
    nm = char(trackResults(t).name);
    rows(end+1,:) = { sprintf('tracking_%s_femur',nm), 'femur', trackResults(t).femur.nFrames, ...
        trackResults(t).femur.posRms_mm, trackResults(t).femur.posP95_mm, trackResults(t).femur.posMax_mm };
    rows(end+1,:) = { sprintf('tracking_%s_tibia',nm), 'tibia', trackResults(t).tibia.nFrames, ...
        trackResults(t).tibia.posRms_mm, trackResults(t).tibia.posP95_mm, trackResults(t).tibia.posMax_mm };
end
rows(end+1,:) = { 'overall','sigma_within_touch',     NaN, sw,    NaN, NaN };
rows(end+1,:) = { 'overall','sigma_between_session',  NaN, sd,    NaN, NaN };
rows(end+1,:) = { 'overall','sigma_tracking',         NaN, sst,   NaN, NaN };
rows(end+1,:) = { 'overall','sigma_combined',         NaN, scomb, NaN, NaN };

Tcsv = cell2table(rows, 'VariableNames', ...
    {'source','landmark','n','rms_mm','p95_mm','max_mm'});
writetable(Tcsv, outCsv);
end

% -------------------------------------------------------------------------
function y = ternary(cond, a, b), if cond, y=a; else, y=b; end, end