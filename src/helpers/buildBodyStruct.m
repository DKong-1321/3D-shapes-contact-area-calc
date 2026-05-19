function body = buildBodyStruct(F, V)
% buildBodyStruct
% ---------------
% Build a struct for a rigid body mesh with:
%   F          : faces (Nx3)
%   V          : vertices (Mx3)
%   triCentroid: per-face centroid
%   triNormal  : per-face unit normal
%   triArea    : per-face area
%   bbox       : 3x2 [xmin xmax; ymin ymax; zmin zmax]
%   kdtree     : KD-tree on triangle centroids
%
% Also creates alias fields:
%   centroids, normals, areas

F = double(F);
V = double(V);

body.F = F;
body.V = V;

% Triangle vertices
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

% Centroids
triCentroid = (v1 + v2 + v3) / 3;
body.triCentroid = triCentroid;
body.centroids   = triCentroid;

% Normals and areas
e1 = v2 - v1;
e2 = v3 - v1;
n  = cross(e1, e2, 2);

area = 0.5 * sqrt(sum(n.^2,2));
body.triArea = area;
body.areas   = area;

nz = vecnorm(n,2,2);
nz(nz == 0) = 1;
triNormal = n ./ nz;

body.triNormal = triNormal;
body.normals   = triNormal;

% Bounding box
xmin = min(V(:,1)); xmax = max(V(:,1));
ymin = min(V(:,2)); ymax = max(V(:,2));
zmin = min(V(:,3)); zmax = max(V(:,3));

body.bbox = [xmin xmax; ymin ymax; zmin zmax];

% -------- ADD THIS (CRITICAL) --------
% KD-tree on triangle centroids for fast neighbour search
body.kdtree = KDTreeSearcher(triCentroid);
% ------------------------------------

end