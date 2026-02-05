% src/computeContactArea_STS_hybrid.m
function results = computeContactArea_STS_hybrid(master, slave, tol, opts)
if nargin < 3 || isempty(tol), tol = 0.2; end
if nargin < 4, opts = struct(); end

master = ensureBodyFields_fast(master);
slave  = ensureBodyFields_fast(slave);

if ~isfield(master,'faceNormals') || isempty(master.faceNormals)
    master.faceNormals = faceNormalsFromFV(master.F, master.V);
end
if ~isfield(master,'kdtreeFaces') || isempty(master.kdtreeFaces)
    master.kdtreeFaces = KDTreeSearcher(master.centroids);
end

[Apen, Atol, maskPen, maskTol] = computeHybridArea(master, slave.F, slave.V, tol);

results = struct();
results.contactArea    = Atol;     % hybrid (penetration + tolerance band)
results.penetrationArea = Apen;    % pure penetration only
results.contactMask    = maskTol;  % mask for hybrid
results.penMask        = maskPen;  % mask for penetration only
results.tolMask        = maskTol;  % same as contactMask here
end


function [Apen, Atol, maskPen, maskTol] = computeHybridArea(master, Fslave, Vslave, tol)
P1 = Vslave(Fslave(:,1),:);
P2 = Vslave(Fslave(:,2),:);
P3 = Vslave(Fslave(:,3),:);
C  = (P1 + P2 + P3) / 3;

idx = knnsearch(master.kdtreeFaces, C, 'K', 1);

Fm  = master.F(idx,:);
V1m = master.V(Fm(:,1),:);
Nm  = master.faceNormals(idx,:);

d1 = dot(P1 - V1m, Nm, 2);
d2 = dot(P2 - V1m, Nm, 2);
d3 = dot(P3 - V1m, Nm, 2);

Aclip0   = zeros(size(Fslave,1),1);
AclipTol = zeros(size(Fslave,1),1);

for i = 1:size(Fslave,1)
    Aclip0(i)   = clippedAreaTriangleHalfSpace(P1(i,:),P2(i,:),P3(i,:), d1(i),d2(i),d3(i));
    AclipTol(i) = clippedAreaTriangleHalfSpace(P1(i,:),P2(i,:),P3(i,:), d1(i)-tol, d2(i)-tol, d3(i)-tol);
end

Apen    = sum(Aclip0);
Atol    = sum(AclipTol);
maskPen = Aclip0   > 0;
maskTol = AclipTol > 0;
end

function N = faceNormalsFromFV(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
N = cross(v2 - v1, v3 - v1, 2);
nrm = sqrt(sum(N.^2,2));
nrm(nrm==0) = 1;
N = N ./ nrm;
end

function A = clippedAreaTriangleHalfSpace(p1,p2,p3, d1,d2,d3)
pts = [p1; p2; p3];
ds  = [d1; d2; d3];

inside = ds <= 0;
if all(~inside), A = 0; return; end
if all(inside),  A = triArea3D(p1,p2,p3); return; end

poly = [];
for k = 1:3
    k2 = mod(k,3) + 1;
    Pk  = pts(k,:);  dk  = ds(k);
    Pk2 = pts(k2,:); dk2 = ds(k2);

    if dk <= 0
        poly = [poly; Pk]; 
    end

    if (dk <= 0 && dk2 > 0) || (dk > 0 && dk2 <= 0)
        t  = dk / (dk - dk2);
        Pi = Pk + t*(Pk2 - Pk);
        poly = [poly; Pi]; 
    end
end

if size(poly,1) < 3, A = 0; return; end

A = 0;
p0 = poly(1,:);
for j = 2:(size(poly,1)-1)
    A = A + triArea3D(p0, poly(j,:), poly(j+1,:));
end
end

function A = triArea3D(a,b,c)
A = 0.5 * norm(cross(b-a, c-a));
end

function body = ensureBodyFields_fast(body)
if ~isfield(body,'centroids')
    if isfield(body,'triCentroid')
        body.centroids = body.triCentroid;
    else
        error('Body missing centroids/triCentroid.');
    end
end
if ~isfield(body,'areas')
    if isfield(body,'triArea')
        body.areas = body.triArea;
    else
        error('Body missing areas/triArea.');
    end
end
if ~isfield(body,'bbox') || isempty(body.bbox)
    V = body.V;
    body.bbox = [min(V(:,1)) max(V(:,1)); min(V(:,2)) max(V(:,2)); min(V(:,3)) max(V(:,3))];
end
end
