function [A_contact, contactMask, tol_eff, d_min] = computeContactArea_adaptiveGap(slave, master, baseTol, opts)
% Adaptive tolerance = min surface gap + baseTol, then call STS.
% slave, master: structs with fields F (faces), V (vertices). kdtree is optional.

    if ~isfield(master, 'kdtree') || isempty(master.kdtree)
        master.kdtree = KDTreeSearcher(master.V);
    end
    if ~isfield(slave, 'kdtree') || isempty(slave.kdtree)
        slave.kdtree = KDTreeSearcher(slave.V);
    end

    F_s = slave.F;
    V_s = slave.V;

    % Triangle centroids on the slave surface
    c = (V_s(F_s(:,1),:) + V_s(F_s(:,2),:) + V_s(F_s(:,3),:)) / 3;  % nTri x 3

    % Nearest-neighbour distances from slave centroids to master surface
    [~, d] = knnsearch(master.kdtree, c);
    d_min  = min(d);                  % minimum 3D gap in this pose

    tol_eff = d_min + baseTol;        % effective tolerance (gap + overshoot)

    % Call existing STS contact algorithm with this effective tolerance
    [A_contact, contactMask] = computeContactArea_STS(slave, master, tol_eff, opts);
end
