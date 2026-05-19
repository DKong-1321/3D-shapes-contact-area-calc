clear; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;
dataDir     = fullfile(projectRoot,'data');
matPath     = fullfile(dataDir,'cylinders.mat');

assert(exist(matPath,'file')==2, 'Missing: %s', matPath);

S = load(matPath);

trajectories = [];
if isfield(S,'trajectories'), trajectories = S.trajectories; end
if isempty(trajectories) && isfield(S,'Trajectories'), trajectories = S.Trajectories; end
if isempty(trajectories)
    fns = fieldnames(S);
    trajectories = S.(fns{1});
end
if istable(trajectories), trajectories = table2array(trajectories); end
if iscell(trajectories), trajectories = [trajectories{:}]; end

tr = trajectories(1);

K = tr.Kinematics;
if isfield(K,'tibiofemoral') && ~isempty(K.tibiofemoral)
    TF = K.tibiofemoral;
elseif isfield(K,'tibiafemoral') && ~isempty(K.tibiafemoral)
    TF = K.tibiafemoral;
else
    error('No tibiofemoral kinematics found in trajectories(1).Kinematics');
end

if istable(TF)
    firstRow = TF(1,:);
    disp('trajectories(1).Kinematics.tibiofemoral first row (table):');
    disp(firstRow);

    vars = TF.Properties.VariableNames;
    for i = 1:numel(vars)
        fprintf('%s = %.4f\n', vars{i}, firstRow{1,i});
    end
else
    firstRow = TF(1,:);
    disp('trajectories(1).Kinematics.tibiofemoral first row (numeric row):');
    disp(firstRow);
end
