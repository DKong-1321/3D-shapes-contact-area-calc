function Cylinder_contact_area_calc
% Cylinder_contact_area_calc
% Approximate contact area between a solid cylinder and a hollow cylinder
% (pipe) from STL files.
% Assumes this file and both STL files are in the SAME folder.
%
% OVERVIEW OF METHOD
% ------------------
% 1. Load the pipe and solid-cylinder meshes from STL files (faces + vertices).
% 2. Treat the pipe as fixed, and move the solid cylinder via a rigid transform:
%       x_new = R * x_old + t
%    (R = rotation matrix, t = translation vector).
% 3. Automatically align the solid cylinder concentrically with the pipe
%    by matching their geometric centres, then apply a user-defined zOffset.
% 4. For each triangle on the SOLID cylinder surface:
%       - Compute its centroid, normal and area.
%       - Find the nearest PIPE vertex to that centroid.
%       - Project the vector (centroid → nearest pipe vertex) onto the
%         triangle normal to get a signed distance.
%       - If that distance is small (within 'tol') and in front of the face
%         (positive sign), treat that triangle as "in contact".
% 5. Sum the areas of all "contact" triangles → approximate contact area.

clc; close all;

% STL filenames for the two shapes
% (moving solid cylinder and fixed hollow cylinder/pipe)
solidFile = 'cylinder 3 it.STL';
pipeFile  = 'cylinder with hole 3 it.STL';

% How far along the Z axis you want to slide the solid cylinder
% relative to a fully concentric position (in STL units, e.g. mm).
%   zOffset = 0      → fully concentric (centres aligned)
%   zOffset < 0      → push the solid cylinder "down" into the pipe
%   zOffset > 0      → pull it "up" out of the pipe
zOffset = 0;         % <<< EDIT THIS TO CONTROL PENETRATION >>>

% Contact tolerance (same units as STL, e.g. mm).
tol = 0.01;

%% Load STL Meshes
% Convert STL files into:
%   F_solid, F_pipe : faces (each row = 3 vertex indices)
%   V_solid, V_pipe : vertices (N x 3 array of [x y z])

[F_solid, V_solid] = loadStlMesh(solidFile, 'solid cylinder');
[F_pipe,  V_pipe ] = loadStlMesh(pipeFile,  'pipe cylinder');

fprintf('Solid cylinder: %d vertices, %d faces\n', size(V_solid,1), size(F_solid,1));
fprintf('Pipe cylinder : %d vertices, %d faces\n', size(V_pipe,1),  size(F_pipe,1));

%% Concentric alignment + rigid transform
% We:
%   1) Compute the geometric centres of both meshes.
%   2) Translate the solid cylinder so its centre matches the pipe's centre.
%   3) Apply an extra zOffset to control penetration along Z.

% Rotation (keep cylinders parallel to each other)
R = eye(3);

% Geometric centres of each mesh
centerSolid = mean(V_solid, 1);   % [cx_solid, cy_solid, cz_solid]
centerPipe  = mean(V_pipe,  1);   % [cx_pipe,  cy_pipe,  cz_pipe ]

% Translation to make them fully concentric, all three coordinates
deltaCenter = centerPipe - centerSolid;

% Now add the user-defined extra shift along Z only
% to control how far into the pipe the solid cylinder goes
deltaCenter(3) = deltaCenter(3) + zOffset;

% Build translation vector t as column
t = deltaCenter(:);   % [tx; ty; tz]

% Apply rigid transform: x_new = R * x_old + t
V_solid_moved = V_solid * R.' + t.';   % N×3

%% Contact area
% Determine which solid-cylinder triangles are in contact with the pipe
% based on centroid–to–pipe proximity along the triangle normal
% and sum their areas

A_contact = contactAreaOnMovingMesh(F_solid, V_solid_moved, V_pipe, tol);

fprintf('\nEstimated contact area: %.6f (STL units^2)\n\n', A_contact);

%% Visualise result
% Plot both meshes to see their relative position and approximate overlap.

figure('Color','w'); hold on; axis equal; grid on;

% Draw pipe (fixed)
trisurf(F_pipe, V_pipe(:,1), V_pipe(:,2), V_pipe(:,3), ...
    'FaceAlpha',0.3, 'EdgeColor','none', 'FaceColor',[0.6 0.6 1]);

% Draw moved solid cylinder
trisurf(F_solid, V_solid_moved(:,1), V_solid_moved(:,2), V_solid_moved(:,3), ...
    'FaceAlpha',0.7, 'EdgeColor','none', 'FaceColor',[1 0.4 0.4]);

xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Solid Cylinder in Pipe – Contact Area ≈ %.4f', A_contact));
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

    % For triangulation objects 
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


function A = contactAreaOnMovingMesh(Fm, Vm, Vfixed, tol)
%contactAreaOnMovingMesh  Approximate contact area on the moving mesh surface.
%
% INPUTS
%   Fm      : moving-mesh faces (nTri x 3 indices into Vm)
%   Vm      : moving-mesh vertices (N x 3), already transformed (R,t applied)
%   Vfixed  : fixed-mesh vertices (M x 3), fixed in space
%   tol     : distance tolerance for deciding "contact" (units of STL)
%
% METHOD
%   For each moving triangle:
%     1. Get its three vertices p1, p2, p3.
%     2. Compute centroid C.
%     3. Compute triangle normal N and area Atri.
%     4. Find nearest fixed-mesh vertex to C.
%     5. Project the vector (C → nearest fixed vertex) onto N to get a
%        signed distance d_signed.
%     6. If |d_signed| < tol AND d_signed > 0, treat this triangle as
%        being in contact with the fixed mesh.
%   Sum the areas of all such triangles.

    % Triangle vertices in 3D for every face
    p1 = Vm(Fm(:,1),:);
    p2 = Vm(Fm(:,2),:);
    p3 = Vm(Fm(:,3),:);

    % Centroids of triangles (average of the 3 vertices)
    C = (p1 + p2 + p3) / 3;

    % Triangle normals
    N = cross(p2 - p1, p3 - p1, 2);

    % Normalise each normal vector to unit length
    normN = vecnorm(N, 2, 2);
    normN(normN == 0) = eps;    % avoid division by zero
    N = N ./ normN;

    % Triangle areas using cross-product magnitude:
    %   area = 0.5 * | (p2-p1) x (p3-p1) |
    Atri = 0.5 * vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);

    % Distances from centroids to all fixed vertices (pairwise, no toolbox)
    % dx, dy, dz are [nTri x nVertFixed] arrays
    dx = C(:,1) - Vfixed(:,1).';
    dy = C(:,2) - Vfixed(:,2).';
    dz = C(:,3) - Vfixed(:,3).';
    D2 = dx.^2 + dy.^2 + dz.^2;   % squared distance

    % For each moving triangle, find index of nearest fixed vertex
    [~, idxMin] = min(D2, [], 2);
    nearestV = Vfixed(idxMin, :);     % coordinates of nearest vertex [nTri x 3]

    % Signed distance along each triangle normal:
    %   v = (centroid → nearest fixed vertex)
    %   d_signed = dot(v, N)
    % Positive d_signed means the fixed vertex lies in front of the face.
    v = nearestV - C;
    d_signed = dot(v, N, 2);

    % Contact triangles:
    %   - close to fixed mesh along normal direction (|d| < tol)
    %   - fixed mesh is in front (d > 0)
    isContact = (abs(d_signed) < tol) & (d_signed > 0);

    % Sum of areas of contacting triangles = approximate contact area
    A = sum(Atri(isContact));

    fprintf('Contact triangles: %d\n', sum(isContact));
end
