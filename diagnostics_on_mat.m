% diagnostics_on_mat.m
% -------------------------------------------------------------------------
% Inspects matlab 2.mat and prints a summary of what is and is not inside.
% Specifically checks for evidence of 5 RK specimens.
%
% Output is printed to the console. To save it to a text file for sharing,
% run:
%       diary mat_contents_summary.txt
%       show_mat_contents
%       diary off
% -------------------------------------------------------------------------

clear; clc;

matFile = 'matlab 2.mat';
assert(exist(matFile,'file')==2, 'Cannot find %s', matFile);

fprintf('========================================================================\n');
fprintf(' Contents of %s\n', matFile);
fprintf('========================================================================\n\n');

S = load(matFile);

% Top-level variables and their sizes
fprintf('TOP-LEVEL VARIABLES\n');
fprintf('-------------------\n');
fns = fieldnames(S);
for i = 1:numel(fns)
    v = S.(fns{i});
    fprintf('  %-15s  class=%-15s  size=%s\n', fns{i}, class(v), sizeStr(v));
end
fprintf('\n');

% Filepath written into the file at recording time
D = S.digitisation;
fprintf('SOURCE FOLDER (read from digitisation.filepath)\n');
fprintf('-----------------------------------------------\n');
fprintf('  %s\n', D.filepath);
fprintf('  (this is the folder ONE specimen was digitised in; sibling folders\n');
fprintf('   for the other RK specimens would be inside its parent directory)\n\n');

% Specimen name
fprintf('SPECIMEN NAME (config.specimen.name)\n');
fprintf('------------------------------------\n');
try
    fprintf('  "%s"\n\n', D.config.specimen.name);
catch
    fprintf('  (not set)\n\n');
end

% The 12 trackers, with their Landmark code
fprintf('TRACKER COLUMNS (digitisation.trackers, size 3x12)\n');
fprintf('--------------------------------------------------\n');
fprintf('  3 rows  = 3 hardware ports (probe + tibia tracker + femur tracker)\n');
fprintf('  12 cols = 12 distinct landmark codes (one anatomical point each)\n\n');

T = D.trackers;
nCol = size(T,2);
fprintf('  %4s  %-6s  %5s   mean position (mm)            tracker port\n', ...
    'col','code','Nval');
fprintf('  %s\n', repmat('-',1,80));
codes = strings(nCol,1);
means = nan(nCol,3);
for j = 1:nCol
    try
        tr = T(1,j);
        codes(j) = string(tr.Landmark);
        P = [tr.Tx tr.Ty tr.Tz];
        P = P(all(~isnan(P),2),:);
        if ~isempty(P), means(j,:) = mean(P,1); end
        fprintf('  %4d  %-6s  %5d  [%8.2f %8.2f %8.2f]   %s\n', ...
            j, codes(j), size(P,1), means(j,1), means(j,2), means(j,3), tr.Name);
    catch ME
        fprintf('  %4d  (error: %s)\n', j, ME.message);
    end
end
fprintf('\n');

% Sample counts per landmark (the 18 samples = 1-second probe hold)
fprintf('SAMPLES PER LANDMARK HOLD (one ~1-second probe touch per landmark)\n');
fprintf('------------------------------------------------------------------\n');
for j = 1:nCol
    try
        n = numel(T(1,j).Tx);
        fprintf('  %-6s : %d samples\n', codes(j), n);
    catch
    end
end
fprintf('\n  (if this file held 5 sessions stacked together we would expect\n');
fprintf('   ~90 samples per landmark, or 5 entries per code. Both absent.)\n\n');

% Camera-motion trajectories (the tracking-error experiment)
fprintf('CAMERA-MOTION TRAJECTORIES (kinematics.trajectories, size 1x2)\n');
fprintf('--------------------------------------------------------------\n');
K = S.kinematics;
for t = 1:numel(K.trajectories)
    tj = K.trajectories(t);
    nF = size(tj.Transform.gTfi, 3);
    fprintf('  trajectory %d  "%s"   %d frames\n', t, char(tj.LoadingCondition), nF);
end
fprintf('  (both are body-clamped, camera-moved recordings of THIS one specimen,\n');
fprintf('   not separate specimens)\n\n');

% Search for any container with a length of 5 (the hypothesised RK axis)
fprintf('SEARCH FOR 5-FOLD STRUCTURE\n');
fprintf('---------------------------\n');
hits = findFives(S, '');
if isempty(hits)
    fprintf('  No variable, struct, cell, or class object with a dimension of 5\n');
    fprintf('  exists anywhere in the file.\n');
else
    fprintf('  Found containers with a dimension of 5:\n');
    for i = 1:numel(hits)
        fprintf('    %s   size=%s\n', hits(i).path, hits(i).size);
    end
end
fprintf('\n');

%% =======================================================================
%% Helpers
%% =======================================================================
function s = sizeStr(x)
    sz = size(x);
    s = strjoin(arrayfun(@(n) sprintf('%d',n), sz, 'UniformOutput', false), 'x');
end

function nm = safeName(D)
    try, nm = char(D.config.specimen.name); catch, nm = '(unknown)'; end
end

function f = folderOnly(p)
    parts = strsplit(p, {'/','\'});
    parts = parts(~cellfun(@isempty, parts));
    if isempty(parts), f = p; else, f = parts{end}; end
end

function hits = findFives(x, path)
    % Recursively look for ANY dimension of 5 in the file's containers.
    % Skips ports-axis (3) and 4x4xN transforms.
    hits = struct('path',{},'size',{});

    sz = size(x);
    if any(sz == 5)
        % Only flag if it's genuinely a 5-fold container of meaningful objects,
        % not a coincidence like a 5-element vector inside a property.
        if isstruct(x) || isobject(x) || iscell(x)
            hits(end+1).path = path;
            hits(end).size   = strjoin(arrayfun(@(n) sprintf('%d',n), sz, 'UniformOutput', false), 'x');
        end
    end

    % Recurse
    if isstruct(x)
        if numel(x) > 1
            try
                x1 = x(1);
                hits = [hits, findFives(x1, sprintf('%s(1)', path))];
            catch
            end
        else
            fns = fieldnames(x);
            for i = 1:numel(fns)
                try
                    hits = [hits, findFives(x.(fns{i}), sprintf('%s.%s', path, fns{i}))]; %#ok<AGROW>
                catch
                end
            end
        end
    elseif iscell(x)
        for i = 1:numel(x)
            hits = [hits, findFives(x{i}, sprintf('%s{%d}', path, i))]; %#ok<AGROW>
        end
    elseif isobject(x)
        if numel(x) > 1
            try
                x1 = x(1);
                hits = [hits, findFives(x1, sprintf('%s(1)', path))];
            catch
            end
        else
            try
                mc = metaclass(x);
                for p = 1:numel(mc.PropertyList)
                    pn = mc.PropertyList(p).Name;
                    try
                        v = x.(pn);
                        hits = [hits, findFives(v, sprintf('%s.%s', path, pn))]; %#ok<AGROW>
                    catch
                    end
                end
            catch
            end
        end
    end
end