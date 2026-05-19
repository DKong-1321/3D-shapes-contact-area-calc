%% debug_alignment_matrixOffsets_showAll.m
% -------------------------------------------------------------------------
% MATRIX-ONLY alignment debug (Jonathan-style)
%
% Uses initial kinematics row to build a constant correction transform:
%   flexion  -> rotation about X
%   anterior -> translation along Y
%
% Then applies it to the 4x4 transform sequence in all plausible ways:
%   - use inv(fTt) or fTt
%   - pre-multiply or post-multiply
%
% Shows ALL cases (no ranking) for frame 1 and a chosen frame.
% -------------------------------------------------------------------------

clear; clc; close all;

%% PATHS
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;
addpath(genpath(fullfile(projectRoot,'src')));

dataDir  = fullfile(projectRoot,'data');
modelDir = fullfile(projectRoot,'model','cylinder');

matPath  = fullfile(dataDir,'cylinders.mat');
femurSTL = fullfile(modelDir,'correct_top.STL');
tibiaSTL = fullfile(modelDir,'correct_bottom 1.STL');

assert(exist(matPath,'file')==2,  'Missing cylinders.mat');
assert(exist(femurSTL,'file')==2, 'Missing femur STL');
assert(exist(tibiaSTL,'file')==2, 'Missing tibia STL');

%% OPTIONS
trialIndex = 1;
frameA = 1;
frameB = 18;  % change if you want
alphaFix = 0.2;
alphaMov = 0.35;

% Jonathan said: "we're happy to offset flexion; would never touch anterior in real world"
% But for debugging, we include both since he suggested it for this dataset.
useFlexOffset    = true;
useAnteriorOffset= true;

% Whether to subtract the initial value (bring it to zero) or add it
% Usually you REMOVE the digitisation offset: so you apply negative of the first-row value.
removeOffset = true;

%% LOAD DATA
S  = load(matPath);
tr = S.trajectories(trialIndex);

fTt = tr.Transform.fTt;
nFrames = size(fTt,3);

%% LOAD STLs
Sf = stlread(femurSTL);
St = stlread(tibiaSTL);

Vm0 = Sf.Points;  Fm = Sf.ConnectivityList; % moving (femur)
Vt0 = St.Points;  Ft = St.ConnectivityList; % fixed (tibia)

%% READ FIRST-ROW KINEMATICS
[flex0_deg, ant0_mm] = getFirstFlexAnt(tr);

fprintf('First row kinematics:\n');
fprintf('  flexion  = %.3f deg\n', flex0_deg);
fprintf('  anterior = %.3f mm\n', ant0_mm);

if removeOffset
    flexUse = -flex0_deg;
    antUse  = -ant0_mm;
else
    flexUse =  flex0_deg;
    antUse  =  ant0_mm;
end

%% BUILD OFFSET TRANSFORM (4x4)
Toffset = eye(4);
if useFlexOffset
    Toffset(1:3,1:3) = rotXdeg(flexUse);
end
if useAnteriorOffset
    Toffset(2,4) = antUse;    % +Y anterior (per your axis convention)
end

disp('Toffset ='); disp(Toffset);

%% DEFINE ALL CASES (show all, no ranking)
% Each case defines how we form the motion transform that moves femur relative to tibia
cases = {
    struct('name','A: Traw = inv(fTt), pre: Toffset*Traw',  'useInv',true,  'mode','pre')
    struct('name','B: Traw = inv(fTt), post: Traw*Toffset', 'useInv',true,  'mode','post')
    struct('name','C: Traw = fTt, pre: Toffset*Traw',       'useInv',false, 'mode','pre')
    struct('name','D: Traw = fTt, post: Traw*Toffset',      'useInv',false, 'mode','post')
};

%% PLOT ALL CASES (two frames per case)
for ci = 1:numel(cases)
    C = cases{ci};

    figure('Color','w','Name',C.name);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ---- frame A ----
    nexttile; hold on; axis equal; grid on; view(3);
    title(sprintf('%s | frame %d', C.name, frameA), 'Interpreter','none');

    patch('Faces',Ft,'Vertices',Vt0,'FaceColor',[0.7 0.7 0.7], ...
          'FaceAlpha',alphaFix,'EdgeColor','none');

    VmA = applyTformToVertices( makeT(ci, frameA, fTt, Toffset, C.useInv, C.mode), Vm0 );
    patch('Faces',Fm,'Vertices',VmA,'FaceColor',[0 0.7 0], ...
          'FaceAlpha',alphaMov,'EdgeColor','none');

    xlabel('X'); ylabel('Y'); zlabel('Z'); camlight; lighting gouraud;

    % ---- frame B ----
    nexttile; hold on; axis equal; grid on; view(3);
    title(sprintf('%s | frame %d', C.name, frameB), 'Interpreter','none');

    patch('Faces',Ft,'Vertices',Vt0,'FaceColor',[0.7 0.7 0.7], ...
          'FaceAlpha',alphaFix,'EdgeColor','none');

    VmB = applyTformToVertices( makeT(ci, frameB, fTt, Toffset, C.useInv, C.mode), Vm0 );
    patch('Faces',Fm,'Vertices',VmB,'FaceColor',[0 0.7 0], ...
          'FaceAlpha',alphaMov,'EdgeColor','none');

    xlabel('X'); ylabel('Y'); zlabel('Z'); camlight; lighting gouraud;
end

fprintf('\nShown %d cases. Pick the one that matches the experiment visually.\n', numel(cases));

%% ====================== HELPERS ======================

function T = makeT(~, k, fTt, Toffset, useInv, mode)
    Traw = fTt(:,:,k);
    if useInv
        Traw = inv(Traw);
    end

    if strcmpi(mode,'pre')
        T = Toffset * Traw;
    else
        T = Traw * Toffset;
    end
end

function [flex0, ant0] = getFirstFlexAnt(tr)
    TF = tr.Kinematics.tibiofemoral;
    if istable(TF)
        vars = TF.Properties.VariableNames;
        flex0 = TF{1, find(strcmpi(vars,'flexion'),1)};
        ant0  = TF{1, find(strcmpi(vars,'anterior'),1)};
    else
        row = TF(1,:);
        flex0 = row(1);
        ant0  = row(5);
    end
end

function Vout = applyTformToVertices(T,V)
    Vh = [V ones(size(V,1),1)];
    Vt = (T*Vh')';
    Vout = Vt(:,1:3);
end

function R = rotXdeg(deg)
    a = deg2rad(deg);
    R = [1 0 0;
         0 cos(a) -sin(a);
         0 sin(a)  cos(a)];
end
