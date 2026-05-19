function varargout = computeContactArea_STS(master, slave, varargin)
% computeContactArea_STS (FAST + ROI prefilter)
% Surface-to-surface static contact area estimate using triangle meshes.
%
% Calling styles:
%   results = computeContactArea_STS(master, slave, opts)
%   [A_contact, contactMask] = computeContactArea_STS(master, slave, tol, opts)
%
% Required fields (from buildBodyStruct):
%   body.F, body.V, body.centroids, body.areas, body.bbox
%
% Optional (will be created if missing):
%   body.kdtree
%
% Speed features:
% - bbox ROI culling (slave centroids)
% - optional additional centroid-distance prefilter (slave centroid must be near master centroid cloud)
% - KNN (knnsearch) + radius filter instead of rangesearch
% - optional centroid-only sampling
% - optional no-mask mode (returnMask=false)

    % ---- Parse inputs ----
    tol  = [];
    opts = struct();

    if numel(varargin) == 1
        if isstruct(varargin{1})
            opts = varargin{1};
        elseif isnumeric(varargin{1}) && isscalar(varargin{1})
            tol = varargin{1};
        else
            error('Invalid third argument. Expected opts struct or scalar tol.');
        end
    elseif numel(varargin) == 2
        if ~(isnumeric(varargin{1}) && isscalar(varargin{1}))
            error('If 4 inputs are used, the 3rd must be scalar tol.');
        end
        tol = varargin{1};
        if isstruct(varargin{2})
            opts = varargin{2};
        else
            error('If 4 inputs are used, the 4th must be opts struct.');
        end
    elseif numel(varargin) > 2
        error('Too many input arguments.');
    end

    % ---- Defaults (chosen for speed) ----
    if ~isfield(opts,'tol'),               opts.tol               = 0.2; end
    if ~isfield(opts,'roiExpandFactor'),   opts.roiExpandFactor   = 1.5; end
    if ~isfield(opts,'neighRadiusFactor'), opts.neighRadiusFactor = 5.0; end
    if ~isfield(opts,'maxNeighbours'),     opts.maxNeighbours     = 25;  end

    % Sampling: 'centroid' (fast) or 'adaptive' (centroid+midpoints if area>threshold)
    if ~isfield(opts,'sampleMode'),        opts.sampleMode        = 'centroid'; end
    if ~isfield(opts,'sampleThreshold'),   opts.sampleThreshold   = 0.5; end  % only used if adaptive

    % Return mask?
    if ~isfield(opts,'returnMask'),        opts.returnMask        = (nargout > 1); end

    % Extra prefilter: keep only candidate slave triangles whose centroids are
    % within (prefilterFactor*tol) of the master centroid cloud. This is very cheap
    % and can massively reduce work for big meshes.
    if ~isfield(opts,'useCentroidPrefilter'), opts.useCentroidPrefilter = true; end
    if ~isfield(opts,'prefilterFactor'),      opts.prefilterFactor      = 6.0;  end

    if ~isempty(tol)
        opts.tol = tol;
    end
    tol = opts.tol;

    % ---- Minimal field checks (fast) ----
    master = ensureBodyFields_fast(master);
    slave  = ensureBodyFields_fast(slave);

    nSlaveTris = size(slave.F, 1);

    % ---- Build KD-tree once if missing ----
    if ~isfield(master,'kdtree') || isempty(master.kdtree)
        master.kdtree = KDTreeSearcher(master.centroids);
    end

    % ---- ROI culling (bbox overlap expanded by tol) ----
    roiMin = max(slave.bbox(:,1), master.bbox(:,1));
    roiMax = min(slave.bbox(:,2), master.bbox(:,2));
    roiMin = roiMin - opts.roiExpandFactor * tol;
    roiMax = roiMax + opts.roiExpandFactor * tol;

    cS = slave.centroids;
    inROI = cS(:,1) >= roiMin(1) & cS(:,1) <= roiMax(1) & ...
            cS(:,2) >= roiMin(2) & cS(:,2) <= roiMax(2) & ...
            cS(:,3) >= roiMin(3) & cS(:,3) <= roiMax(3);

    candidateIdx = find(inROI);
    if isempty(candidateIdx)
        A_contact = 0;
        if opts.returnMask
            contactMask = false(nSlaveTris,1);
        else
            contactMask = [];
        end
        [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
        return
    end

    % ---- Additional centroid-distance prefilter (optional) ----
    neighRadius = max(1e-12, opts.neighRadiusFactor * tol);
    if opts.useCentroidPrefilter
        % nearest master centroid distance for each candidate slave centroid
        cCand = cS(candidateIdx,:);
        % K=1 is enough for prefilter
        [~, d0] = knnsearch(master.kdtree, cCand, 'K', 1);

        % allow a bit more than neighRadius to avoid false negatives
        gate = max(neighRadius, opts.prefilterFactor * tol);
        keep = (d0 <= gate);

        candidateIdx = candidateIdx(keep);

        if isempty(candidateIdx)
            A_contact = 0;
            if opts.returnMask
                contactMask = false(nSlaveTris,1);
            else
                contactMask = [];
            end
            [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
            return
        end
    end

    % ---- Precompute master triangle vertex arrays (big speedup) ----
    Fm = master.F; Vm = master.V;
    Am = Vm(Fm(:,1),:);
    Bm = Vm(Fm(:,2),:);
    Cm = Vm(Fm(:,3),:);

    Fs = slave.F;  Vs = slave.V;

    if opts.returnMask
        contactMask = false(nSlaveTris, 1);
    else
        contactMask = [];
    end

    K = max(1, opts.maxNeighbours);

    % ---- Main loop ----
    A_contact = 0;

    for ii = 1:numel(candidateIdx)
        tIdx = candidateIdx(ii);

        % pick sample points (default: centroid only)
        sTri = Fs(tIdx,:);
        p1 = Vs(sTri(1),:);
        p2 = Vs(sTri(2),:);
        p3 = Vs(sTri(3),:);

        ps = samplePointsFast(p1,p2,p3, slave.areas(tIdx), opts);

        hit = false;
        for sp = 1:size(ps,1)
            q = ps(sp,:);

            % KNN then radius filter (faster than rangesearch)
            [neighIdx, neighDist] = knnsearch(master.kdtree, q, 'K', K);

            inR = neighDist <= neighRadius;
            if ~any(inR)
                continue
            end
            neighIdx = neighIdx(inR);

            % check distance to each neighbour master triangle
            for k = 1:numel(neighIdx)
                j = neighIdx(k);
                d = pointTriangleDistance_fast(q, Am(j,:), Bm(j,:), Cm(j,:));
                if d <= tol
                    hit = true;
                    break
                end
            end

            if hit
                break
            end
        end

        if hit
            A_contact = A_contact + slave.areas(tIdx);
            if opts.returnMask
                contactMask(tIdx) = true;
            end
        end
    end

    [varargout{1:nargout}] = packOutputs(A_contact, contactMask);
end

% ===================== Helpers =====================

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

function ps = samplePointsFast(p1,p2,p3, triArea, opts)
    c  = (p1+p2+p3)/3;

    if strcmpi(opts.sampleMode,'adaptive') && triArea > opts.sampleThreshold
        m12 = (p1+p2)/2;
        m23 = (p2+p3)/2;
        m31 = (p3+p1)/2;
        ps = [c; m12; m23; m31];
    else
        ps = c;
    end
end

function d = pointTriangleDistance_fast(p, a, b, c)
% Returns Euclidean distance from point p to triangle abc in 3D.
    ab = b - a;
    ac = c - a;
    ap = p - a;

    d1 = dot(ab, ap);
    d2 = dot(ac, ap);
    if d1 <= 0 && d2 <= 0
        d = norm(ap); return;
    end

    bp = p - b;
    d3 = dot(ab, bp);
    d4 = dot(ac, bp);
    if d3 >= 0 && d4 <= d3
        d = norm(bp); return;
    end

    vc = d1*d4 - d3*d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3);
        proj = a + v * ab;
        d = norm(p - proj); return;
    end

    cp = p - c;
    d5 = dot(ab, cp);
    d6 = dot(ac, cp);
    if d6 >= 0 && d5 <= d6
        d = norm(cp); return;
    end

    vb = d5*d2 - d1*d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6);
        proj = a + w * ac;
        d = norm(p - proj); return;
    end

    va = d3*d6 - d5*d4;
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        proj = b + w * (c - b);
        d = norm(p - proj); return;
    end

    denom = 1 / (va + vb + vc);
    v = vb * denom;
    w = vc * denom;
    proj = a + ab * v + ac * w;
    d = norm(p - proj);
end

function varargout = packOutputs(A_contact, contactMask)
    if nargout <= 1
        results = struct();
        results.contactArea = A_contact;
        results.contactMask = contactMask;
        varargout = {results};
    else
        varargout = {A_contact, contactMask};
    end
end