function Cube_pyramid_area_calc
% Cube_pyramid_area_calc
% Approximate contact area between a pyramid and a cube from STL files.
% Assumes this file and both STL files are in the SAME folder.
%
% OVERVIEW OF METHOD
% ------------------
% 1. Load the cube and pyramid meshes from STL files (faces + vertices).
% 2. Treat the cube as fixed, and move the pyramid via a rigid transform:
%       x_new = R * x_old + t
%    (R = rotation matrix, t = translation vector).
% 3. For each triangle on the PYRAMID surface:
%       - Compute its centroid and area.
%       - Compute the shortest distance from the centroid to the cube's
%         axis-aligned bounding box.
%       - If that distance is small (within 'tol'), treat that triangle
%         as "in contact".
% 4. Sum the areas of all "contact" triangles → approximate contact area.
% 5. Visualise contact vs non-contact triangles in different colours.

clc; close all;

% STL filenames for the two shapes
pyrFile  = '50mm height pyramid 4 iterations.stl'; % moving object pyramid
cubeFile = '50mm cube 4 iterations.stl';           % fixed object cube

% Rigid transform applied to the PYRAMID only:
%   - R is a 3x3 rotation matrix
%   - t is a 3x1 translation vector [tx; ty; tz]
% Here: no rotation, just move the pyramid down by 5 units in Z.
R   = eye(3);      % 3x3 rotation identity matrix (no rotation)
t   = [0; 0; -5];  % translation vector [x; y; z] - move pyramid into cube 

% Contact tolerance (same units as STL, e.g. mm).
% Any pyramid triangle whose centroid is within 'tol' of the cube's
% bounding box is counted as "in contact".
tol = 0.5;         % try 0.5–1.0 mm; tune as needed

%% Load STL Mesh
% Convert STL files into:
%   Fp, Fc : faces (each row = 3 vertex indices)
%   Vp, Vc : vertices (N x 3 array of [x y z])

[Fp, Vp] = loadStlMesh(pyrFile,  'pyramid');
[Fc, Vc] = loadStlMesh(cubeFile, 'cube');

fprintf('Pyramid: %d vertices, %d faces\n', size(Vp,1), size(Fp,1));
fprintf('Cube   : %d vertices, %d faces\n', size(Vc,1), size(Fc,1));

%% Rigid Transform of pyramid
% Rigid body transform:
%   x_new = R * x_old + t
% Implemented in row-wise form as:
%   Vp_moved = Vp * R.' + t.'

Vp_moved = Vp * R.' + t.';  

%% Contact area and mask
[A_contact, isContact] = contactAreaOnPyramid(Fp, Vp_moved, Vc, tol);

fprintf('\nEstimated contact area: %.6f (STL units^2)\n\n', A_contact);

%% Result visualisation
figure('Color','w'); hold on; axis equal; grid on;

% Draw cube (fixed)
trisurf(Fc, Vc(:,1), Vc(:,2), Vc(:,3), ...
    'FaceAlpha',0.3, 'EdgeColor','none', 'FaceColor',[0.6 0.6 1]);

% Split pyramid faces into contact / non-contact
contactFaces    = Fp(isContact, :);
nonContactFaces = Fp(~isContact, :);

% Draw non-contact pyramid faces (light red)
if ~isempty(nonContactFaces)
    trisurf(nonContactFaces, ...
        Vp_moved(:,1), Vp_moved(:,2), Vp_moved(:,3), ...
        'FaceAlpha',0.7, 'EdgeColor','none', 'FaceColor',[1.0 0.6 0.6]);
end

% Draw contact pyramid faces on top (green)
if ~isempty(contactFaces)
    trisurf(contactFaces, ...
        Vp_moved(:,1), Vp_moved(:,2), Vp_moved(:,3), ...
        'FaceAlpha',1.0, 'EdgeColor','none', 'FaceColor',[0.2 0.8 0.2]);
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Pyramid into Cube – Contact Area ≈ %.4f', A_contact));
view(3);
legend({'Cube','Pyramid (non-contact)','Pyramid (contact)'}, 'Location','bestoutside');

end  % end main function


%% --------- Helper functions --------- %

function [F, V] = loadStlMesh(fname, nameStr)
%loadStlMesh  Load STL as faces F and vertices V (N×3) using triangulation.

    TR = stlread(fname);   % should return a triangulation or similar object

    % For triangulation objects (MATLAB's built-in stlread)
    if isa(TR, 'triangulation')
        F = TR.ConnectivityList;
        V = TR.Points;
    else
        % If your stlread returns a struct with fields, handle that here:
        if isfield(TR, 'ConnectivityList') && isfield(TR, 'Points')
            F = TR.ConnectivityList;
            V = TR.Points;
        elseif isfield(TR, 'faces') && isfield(TR, 'vertices')
            F = TR.faces;
            V = TR.vertices;
        else
            error('Unsupported STL format returned by stlread for %s.', nameStr);
        end
    end

    % Make sure V is N×3
    V = ensureNx3(V, nameStr);
end


function V = ensureNx3(Vraw, nameStr)
%ensureNx3  Make sure vertex array is N×3.

    Vraw = double(Vraw);
    [nr, nc] = size(Vraw);

    if nc == 3
        V = Vraw;

    elseif nr == 3
        V = Vraw.';

    else
        ne = numel(Vraw);
        if mod(ne, 3) ~= 0
            error('%s vertices have unexpected size %d×%d (total %d not divisible by 3).', ...
                  nameStr, nr, nc, ne);
        end
        V = reshape(Vraw, [], 3);
        fprintf('Reshaped %s vertices from %d×%d to %d×3\n', ...
                nameStr, nr, nc, size(V,1));
    end
end


function [A, isContact] = contactAreaOnPyramid(Fp, Vp, Vc, tol)
%contactAreaOnPyramid  Approximate contact area on the pyramid surface.
%
% INPUTS
%   Fp  : pyramid faces (nTri x 3 indices into Vp)
%   Vp  : pyramid vertices (N x 3), already transformed (R,t applied)
%   Vc  : cube vertices (M x 3), fixed in space
%   tol : distance tolerance for deciding "contact" (units of STL)
%
% OUTPUTS
%   A         : total contact area (sum of areas of contact triangles)
%   isContact : logical vector (nTri x 1), true where face is in contact
%
% METHOD
%   For each pyramid triangle:
%     1. Compute centroid C and triangle area Atri.
%     2. Compute closest point Q on cube's axis-aligned bounding box.
%     3. Distance d = ||C - Q||.
%     4. If d <= tol, mark as contact.

    % Triangle vertices
    p1 = Vp(Fp(:,1),:);
    p2 = Vp(Fp(:,2),:);
    p3 = Vp(Fp(:,3),:);

    % Centroids
    C = (p1 + p2 + p3) / 3;      % [nTri x 3]
    cx = C(:,1); cy = C(:,2); cz = C(:,3);

    % Triangle areas
    Atri = 0.5 * vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);

    % Cube bounding box
    xmin = min(Vc(:,1)); xmax = max(Vc(:,1));
    ymin = min(Vc(:,2)); ymax = max(Vc(:,2));
    zmin = min(Vc(:,3)); zmax = max(Vc(:,3));

    % Closest point on the axis-aligned box to each centroid
    qx = min(max(cx, xmin), xmax);
    qy = min(max(cy, ymin), ymax);
    qz = min(max(cz, zmin), zmax);

    % Distances from centroids to the box
    dx = cx - qx;
    dy = cy - qy;
    dz = cz - qz;
    d  = sqrt(dx.^2 + dy.^2 + dz.^2);

    % Contact mask: within tol of the cube's bounding box
    isContact = (d <= tol);

    % Contact area
    A = sum(Atri(isContact));

    fprintf('Contact triangles: %d (out of %d)\n', sum(isContact), numel(isContact));
end
