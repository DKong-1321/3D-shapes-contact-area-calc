function results = computeContactArea_STS_penetration(master, slave, varargin)

% ---------------- INPUT PARSE ----------------
tol  = [];
opts = struct();

if nargin == 3
    if isstruct(varargin{1})
        opts = varargin{1};
    else
        tol = varargin{1};
    end
elseif nargin == 4
    tol  = varargin{1};
    opts = varargin{2};
end

% ---------------- DEFAULTS ----------------
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

% ---------------- ENSURE FIELDS ----------------
master = ensureBodyFields_fast(master);
slave  = ensureBodyFields_fast(slave);

nSlaveTris = size(slave.F,1);

if ~isfield(master,'kdtree') || isempty(master.kdtree)
    master.kdtree = KDTreeSearcher(master.centroids);
end

% ---------------- ROI ----------------
roiMin = max(slave.bbox(:,1), master.bbox(:,1)) - opts.roiExpandFactor*tol;
roiMax = min(slave.bbox(:,2), master.bbox(:,2)) + opts.roiExpandFactor*tol;

cS = slave.centroids;
inROI = all(cS >= roiMin' & cS <= roiMax',2);
candidateIdx = find(inROI);

if isempty(candidateIdx)
    results = emptyResults(nSlaveTris); return
end

% ---------------- PREFILTER ----------------
neighRadius = max(1e-12, opts.neighRadiusFactor*tol);
d0 = [];

if opts.useCentroidPrefilter
    [~, d0] = knnsearch(master.kdtree, cS(candidateIdx,:), 'K', 1);
    keep = d0 <= max(neighRadius, opts.prefilterFactor*tol);
    candidateIdx = candidateIdx(keep);
    d0 = d0(keep);
    if isempty(candidateIdx)
        results = emptyResults(nSlaveTris); return
    end
end

% ---------------- PRECOMPUTE ----------------
Fm = master.F; Vm = master.V;
Am = Vm(Fm(:,1),:); Bm = Vm(Fm(:,2),:); Cm = Vm(Fm(:,3),:);
Fs = slave.F;  Vs = slave.V;

contactMask = false(nSlaveTris,1);
penHit = false(numel(candidateIdx),1);

% ---------------- PENETRATION PASS ----------------
if opts.includePenetration && opts.useInpolyhedron && exist('inpolyhedron','file')==2
    if isempty(d0)
        [~, d0] = knnsearch(master.kdtree, cS(candidateIdx,:), 'K', 1);
    end
    idxNear = d0 <= opts.penetrationGate_mm;
    if any(idxNear)
        candNear = candidateIdx(idxNear);
        [ptsXYZ, triMap] = samplePointsForTriangles( ...
            Fs(candNear,:), Vs, slave.areas(candNear), opts);

        inside = inpolyhedron(Fm, Vm, ptsXYZ);
        penHit(idxNear) = accumarray(triMap, inside, ...
            [numel(candNear) 1], @any, false);
    end
end

% ---------------- MAIN LOOP ----------------
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
