% main.m
% Entry point for the code
clear;
clc;

% Add paths to helpers, transforms, and models folder
addpath(genpath('src/helpers'));
addpath(genpath('src/transforms'));
addpath('model');

% Set file paths for models
cubeFile = 'model/50mm cube 4 iterations.stl';
cylinderFile = 'model/cylinder 30mm diameter 4 it.stl';

% Contact tolerance and options for the calculation
tol = 0.25;  % Define contact tolerance
opts.sampleThreshold   = 0.2;    % fraction of samples within tol to call triangle "in contact"
opts.neighRadiusFactor = 5.0;    % neighbour search radius = factor * tol
opts.maxNeighbours     = 30;     % max number of candidate master tris per sample point
opts.roiExpandFactor   = 1.5;    % expand ROI box by this * tol

% Call the main function to calculate contact area
ContactArea_STL_STS_cube_cylinder();
