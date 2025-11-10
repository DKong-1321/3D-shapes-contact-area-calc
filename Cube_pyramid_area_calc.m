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
%       - Compute its centroid, normal and area.
%       - Find the nearest CUBE vertex to that centroid.
%       - Project the vector (centroid → nearest cube vertex) onto the
%         triangle normal to get a signed distance.
%       - If that distance is small (within 'tol') and in front of the face
%         (positive sign), treat that triangle as "in contact".
% 4. Sum the areas of all "contact" triangles → approximate contact area.

clc; close all;

% STL filenames for the two shapes
pyrFile  = '50mm height pyramid 4 iterations.stl';
cubeFile = '50mm cube 4 iterations.stl';

% Rigid transform applied to the PYRAMID only:
%   - R is a 3x3 rotation matrix
%   - t is a 3x1 translation vector [tx; ty; tz]
% Here: no rotation, just move the pyramid down by 5 units in Z.
R   = eye(3);      % 3x3 rotation matrix (no rotation)
t   = [0; 0; -5];  % translation [x; y; z] - move pyramid down into cube

% Contact tolerance (same units as STL, e.g. mm).
% Any pyramid triangle whose centroid is within 'tol' of the cube (along
% its normal direction) is counted as "in contact".
tol = 0.01;

%% Load STL Mesh
% Convert STL files into:
%   Fp, Fc : faces (each row = 3 vertex indices)
%   Vp, Vc : vertices (N x 3 array of [x y z])

[Fp, Vp] = loadStlMesh(pyrFile,  'pyramid');
[Fc, Vc] = loadStlMesh(cubeFile, 'cube');

fprintf('Pyramid: %d vertices, %d faces\n', size(Vp,1), size(Fp,1));
fprintf('Cube   : %d vertices, %d faces\n', size(Vc,1), size(Fc,1));

%% Rigid Transform of pyramid
% Vp is N×3, R is 3×3, t is 3×1
%
% Rigid body transform:
%   x_new = R * x_old + t
% Implemented in row-wise form as:
%   Vp_moved = Vp * R.' + t.'
%
% Note: the cube vertices Vc are NOT moved. The cube is treated as fixed
% in space; the pyramid is moved relative to it.

Vp_moved = Vp * R.' + t.';   % N×3

%% Contact area
% Determine which pyramid triangles are "in contact" with the cube
% (based on centroid–to–cube proximity along the triangle normal)
% and sum their areas.

A_contact = contactAreaOnPyramid(Fp, Vp_moved, Vc, tol);

fprintf('\nEstimated contact area: %.6f (STL units^2)\n\n', A_contact);

%% Result visualisation
% Plot both meshes to see their relative position and approximate overlap.

figure('Color','w'); hold on; axis equal; grid on;

% Draw cube (fixed)
trisurf(Fc, Vc(:,1), Vc(:,2), Vc(:,3), ...
    'FaceAlpha',0.3, 'EdgeColor','none', 'FaceColor',[0.6 0.6 1]);

% Draw moved pyramid
trisurf(Fp, Vp_moved(:,1), Vp_moved(:,2), Vp_moved(:,3), ...
    'FaceAlpha',0.7, 'EdgeColor','none', 'FaceColor',[1 0.4 0.4]);

xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Pyramid into Cube – Contact Area ≈ %.4f', A_contact));
view(3);

end  % end of main function


%% Functions

function [F, V] = loadStlMesh(fname, nameStr)
%loadStlMesh  Load STL as faces F and vertices V (N×3) using triangulation.
%
% This wraps whatever 'stlread' you have, and normalises its output to:
%   F : nFaces x 3  (indices into V)
%   V : nVerts x 3  (x, y, z coordinates)
%
% We prefer the triangulation-style output (ConnectivityList, Points) to
% avoid weird reshaping issues.

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

    % Make sure V is N×3 (and not 3×N or flattened 1×(3N) etc.)
    V = ensureNx3(V, nameStr);
end


function V = ensureNx3(Vraw, nameStr)
%ensureNx3  Make sure vertex array is N×3.
%
% Accepts:
%   - N×3    → unchanged
%   - 3×N    → transposed to N×3
%   - 1×(3N) or (3N)×1 → reshaped to N×3
%
% The goal is simply to standardise the vertex array shape so the rest of
% the code can assume rows correspond to [x y z].

    Vraw = double(Vraw);       % make sure it's numeric
    [nr, nc] = size(Vraw);

    if nc == 3
        % Already N×3
        V = Vraw;

    elseif nr == 3
        % Was 3×N → make N×3
        V = Vraw.';

    else
        % Flatten and reshape if total elements is a multiple of 3
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


function A = contactAreaOnPyramid(Fp, Vp, Vc, tol)
%contactAreaOnPyramid  Approximate contact area on the pyramid surface.
%
% INPUTS
%   Fp  : pyramid faces (nTri x 3 indices into Vp)
%   Vp  : pyramid vertices (N x 3), already transformed (R,t applied)
%   Vc  : cube vertices (M x 3), fixed in space
%   tol : distance tolerance for deciding "contact" (units of STL)
%
% METHOD
%   For each pyramid triangle:
%     1. Get its three vertices p1, p2, p3.
%     2. Compute centroid C.
%     3. Compute triangle normal N and area Atri.
%     4. Find nearest cube vertex to C.
%     5. Project the vector (C → nearest cube vertex) onto N to get a
%        signed distance d_signed.
%     6. If |d_signed| < tol AND d_signed > 0, treat this triangle as
%        being in contact with the cube.
%   Sum the areas of all such triangles.

    % Triangle vertices in 3D for every face
    p1 = Vp(Fp(:,1),:);
    p2 = Vp(Fp(:,2),:);
    p3 = Vp(Fp(:,3),:);

    % Centroids of triangles (average of the 3 vertices)
    C = (p1 + p2 + p3) / 3;

    % Triangle normals
    N = cross(p2 - p1, p3 - p1, 2);
    % Normalise each normal vector to unit length
    N = N ./ vecnorm(N, 2, 2);

    % Triangle areas using cross-product magnitude:
    %   area = 0.5 * | (p2-p1) x (p3-p1) |
    Atri = 0.5 * vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);

    % Distances from centroids to all cube vertices (pairwise, no toolbox)
    % dx, dy, dz are [nTri x nVertC] arrays
    dx = C(:,1) - Vc(:,1).';
    dy = C(:,2) - Vc(:,2).';
    dz = C(:,3) - Vc(:,3).';
    D2 = dx.^2 + dy.^2 + dz.^2;   % squared distance

    % For each pyramid triangle, find index of nearest cube vertex
    [~, idxMin] = min(D2, [], 2);
    nearestV = Vc(idxMin, :);     % coordinates of nearest vertex [nTri x 3]

    % Signed distance along each triangle normal:
    %   v = (centroid → nearest cube vertex)
    %   d_signed = dot(v, N)
    % Positive d_signed means the cube vertex lies in front of the face.
    v = nearestV - C;
    d_signed = dot(v, N, 2);

    % "Contact" triangles:
    %   - close to cube along normal direction (|d| < tol)
    %   - cube is in front (d > 0)
    isContact = (abs(d_signed) < tol) & (d_signed > 0);

    % Sum of areas of contacting triangles = approximate contact area
    A = sum(Atri(isContact));

    fprintf('Contact triangles: %d\n', sum(isContact));
end
