function results = computeContactArea_STS_penetration(master, slave, varargin)

% Parse inputs
tol  = [];
opts = struct();

if numel(varargin) == 1
    if isstruct(varargin{1})
        opts = varargin{1};
    else
        tol = varargin{1};
    end
elseif numel(varargin) == 2
    tol  = varargin{1};
    opts = varargin{2};
end

% Defaults
if ~isfield(opts,'tol'),                 opts.tol = 0.2; end
if ~isfield(opts,'roiExpandFactor'),     opts.roiExpandFactor = 1.5; end
if ~isfield(opts,'neighRadiusFactor'),   opts.neighRadiusFactor = 5.0; end
if ~isfield(opts,'maxNeighbours'),       opts.maxNeighbours = 25; end
if ~isfield(opts,'sampleMode'),          opts.sampleMode = 'centroid'; end
if ~isfield(opts,'sampleThreshold'),     opts.sampleThreshold = 0.5; end
if ~isfield(opts,'useCentroidPrefilter'),opts.useCentroidPrefilter = true; end
if ~isfield(opts,'prefilterFactor'),     opts.prefilterFactor = 6.0; end
if ~isfield(opts,'includePenetration'),  opts.includePenetration = false; end
if ~isfield(opts,'penetrationGate_mm'),  opts.penetrationGate_mm = 3.0; end
if ~isfield(opts,'useInpolyhedron'),     opts.useInpolyhedron = true; end

if ~isempty(tol), opts.tol = tol; end
tol = opts.tol;

% Ensure fields
master = ensureBodyFields_fast(master);
slave  = ensureBodyFields_fast(slave);

nSlaveTris = size(slave.F,1);

if ~isfield(master,'kdtree') || isempty(master.kdtree)
    master.kdtree = KDTreeSearcher(master.centroids);
end

% ROI
roiMin = max(slave.bbox(:,1), master.bbox(:,1)) - opts.roiExpandFactor*tol;
roiMax = min(slave.bbox(:,2), master.bbox(:,2)) + opts.roiExpandFactor*tol;

cS = slave.centroids;
inROI = cS(:,1)>=roiMin(1) & cS(:,1)<=roiMax(1) & ...
        cS(:,2)>=roiMin(2) & cS(:,2)<=roiMax(2) & ...
        cS(:,3)>=roiMin(3) & cS(:,3)<=roiMax(3);

candidateIdx = find(inROI);
if isempty(candidateIdx)
    results = emptyResults(nSlaveTris); return
end

% Prefilter
neighRadius = max(1e-12, opts.neighRadiusFactor*tol);
d0 = [];

if opts.useCentroidPrefilter
    cCand = cS(candidateIdx,:);
    [~, d0] = knnsearch(master.kdtree, cCand, 'K', 1);
    keep = d0 <= max(neighRadius, opts.prefilterFactor*tol);
    candidateIdx = candidateIdx(keep);
    d0 = d0(keep);
    if isempty(candidateIdx)
        results = emptyResults(nSlaveTris); return
    end
end

% Precompute
Fm = master.F; Vm = master.V;
Am = Vm(Fm(:,1),:); Bm = Vm(Fm(:,2),:); Cm = Vm(Fm(:,3),:);
Fs = slave.F;  Vs = slave.V;

contactMask = false(nSlaveTris,1);
penHit = false(numel(candidateIdx),1);

% Penetration pass
if opts.includePenetration && opts.useInpolyhedron && exist('inpolyhedron','file')==2
    if isempty(d0)
        [~, d0] = knnsearch(master.kdtree, cS(candidateIdx,:), 'K', 1);
    end
    idxNear = d0 <= opts.penetrationGate_mm;
    if any(idxNear)
        candNear = candidateIdx(idxNear);
        [ptsXYZ, triMap] = samplePointsForTriangles(Fs(candNear,:), Vs, slave.areas(candNear), opts);
        inside = inpolyhedron(Fm, Vm, ptsXYZ);
        penHit(idxNear) = accumarray(triMap, inside, [numel(candNear) 1], @any, false);
    end
end

% ---- Main loop ----
A_surface = 0;
A_pen = 0;
K = max(1, opts.maxNeighbours);

for ii = 1:numel(candidateIdx)
    tIdx = candidateIdx(ii);

    if opts.includePenetration && penHit(ii)
        contactMask(tIdx) = true;
        A_pen = A_pen + slave.areas(tIdx);
        continue
    end

    tri = Fs(tIdx,:);
    p1 = Vs(tri(1),:); p2 = Vs(tri(2),:); p3 = Vs(tri(3),:);
    ps = samplePointsFast(p1,p2,p3, slave.areas(tIdx), opts);

    hit = false;
    for sp = 1:size(ps,1)
        q = ps(sp,:);
        [neighIdx, neighDist] = knnsearch(master.kdtree, q, 'K', K);
        neighIdx = neighIdx(neighDist <= neighRadius);
        for k = 1:numel(neighIdx)
            j = neighIdx(k);
            if pointTriangleDistance_fast(q, Am(j,:), Bm(j,:), Cm(j,:)) <= tol
                hit = true; break
            end
        end
        if hit, break; end
    end

    if hit
        contactMask(tIdx) = true;
        A_surface = A_surface + slave.areas(tIdx);
    end
end

results.contactArea = A_surface + A_pen;
results.contactMask = contactMask;
results.contactArea_surfaceOnly = A_surface;
results.contactArea_penetrationOnly = A_pen;

end

% HELPERS

function results = emptyResults(n)
results.contactArea = 0;
results.contactMask = false(n,1);
results.contactArea_surfaceOnly = 0;
results.contactArea_penetrationOnly = 0;
end

function body = ensureBodyFields_fast(body)
if ~isfield(body,'centroids')
    body.centroids = body.triCentroid;
end
if ~isfield(body,'areas')
    body.areas = body.triArea;
end
if ~isfield(body,'bbox') || isempty(body.bbox)
    V = body.V;
    body.bbox = [min(V(:,1)) max(V(:,1)); min(V(:,2)) max(V(:,2)); min(V(:,3)) max(V(:,3))];
end
end

function ps = samplePointsFast(p1,p2,p3, triArea, opts)
c = (p1+p2+p3)/3;
if strcmpi(opts.sampleMode,'adaptive') && triArea > opts.sampleThreshold
    ps = [c; (p1+p2)/2; (p2+p3)/2; (p3+p1)/2];
else
    ps = c;
end
end

function [ptsXYZ, triMap] = samplePointsForTriangles(Ftri, V, triAreas, opts)
nTri = size(Ftri,1);
ptsXYZ = [];
triMap = [];
for i = 1:nTri
    a = V(Ftri(i,1),:); b = V(Ftri(i,2),:); c = V(Ftri(i,3),:);
    cen = (a+b+c)/3;
    if strcmpi(opts.sampleMode,'adaptive') && triAreas(i) > opts.sampleThreshold
        P = [cen; (a+b)/2; (b+c)/2; (c+a)/2];
    else
        P = cen;
    end
    ptsXYZ = [ptsXYZ; P];
    triMap = [triMap; i*ones(size(P,1),1)]; 
end
end

function d = pointTriangleDistance_fast(p,a,b,c)
ab=b-a; ac=c-a; ap=p-a;
d1=dot(ab,ap); d2=dot(ac,ap);
if d1<=0 && d2<=0, d=norm(ap); return; end
bp=p-b; d3=dot(ab,bp); d4=dot(ac,bp);
if d3>=0 && d4<=d3, d=norm(bp); return; end
vc=d1*d4-d3*d2;
if vc<=0 && d1>=0 && d3<=0
    v=d1/(d1-d3); d=norm(p-(a+v*ab)); return
end
cp=p-c; d5=dot(ab,cp); d6=dot(ac,cp);
if d6>=0 && d5<=d6, d=norm(cp); return; end
vb=d5*d2-d1*d6;
if vb<=0 && d2>=0 && d6<=0
    w=d2/(d2-d6); d=norm(p-(a+w*ac)); return
end
va=d3*d6-d5*d4;
if va<=0 && (d4-d3)>=0 && (d5-d6)>=0
    w=(d4-d3)/((d4-d3)+(d5-d6)); d=norm(p-(b+w*(c-b))); return
end
den=1/(va+vb+vc);
v=vb*den; w=vc*den;
d=norm(p-(a+ab*v+ac*w));
end
