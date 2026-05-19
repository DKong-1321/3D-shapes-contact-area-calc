% view_cylinder_stls.m
% -----------------------------------------
% Reads two STL files and visualises them.
% Place in project root and run.
% -----------------------------------------

clear; clc; close all;

% --- Paths (edit if your filenames differ) ---
projectRoot = fileparts(mfilename('fullpath'));
femurStl = fullfile(projectRoot,'model','cylinder','top.stl');
tibiaStl = fullfile(projectRoot,'model','cylinder','bottom.stl');

assert(exist(femurStl,'file')==2, 'Missing: %s', femurStl);
assert(exist(tibiaStl,'file')==2, 'Missing: %s', tibiaStl);

% --- Read STLs ---
Sf = stlread(femurStl);
St = stlread(tibiaStl);

Vf = double(Sf.Points);   Ff = Sf.ConnectivityList;
Vt = double(St.Points);   Ft = St.ConnectivityList;

% --- Optional: centre each STL at its centroid (toggle) ---
doCenter = false;
if doCenter
    Vf = Vf - mean(Vf,1);
    Vt = Vt - mean(Vt,1);
end

% --- Plot ---
figure('Color','w','Name','Cylinder STLs');
ax = axes; hold(ax,'on'); axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'X (mm)'); ylabel(ax,'Y (mm)'); zlabel(ax,'Z (mm)');

hTibia = patch('Faces',Ft,'Vertices',Vt, ...
    'FaceColor',[0.27 0.51 0.71],'FaceAlpha',0.5,'EdgeColor','none');
hFemur = patch('Faces',Ff,'Vertices',Vf, ...
    'FaceColor',[0.94 0.50 0.31],'FaceAlpha',0.9,'EdgeColor','none');

camlight(ax,'headlight');
lighting(ax,'gouraud');

legend([hTibia hFemur], {'Tibia STL','Femur STL'}, 'Location','best');
title(ax, sprintf('Femur faces: %d | Tibia faces: %d', size(Ff,1), size(Ft,1)));
