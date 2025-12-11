function [A_contact, contactMask, tol_eff, e_min] = ...
    computeContactArea_adaptiveGapVertical(slave, master, baseTol, opts)
% Adaptive tolerance = min vertical gap + baseTol, then call STS
% Vertical gap is measured along the global Z axis
% slave, master: structs with fields F (faces), V (vertices)

    % Ensure KD-tree for STS (3D) exists
    if ~isfield(master, 'kdtree') || isempty(master.kdtree)
        master.kdtree = KDTreeSearcher(master.V);
    end
    if ~isfield(slave, 'kdtree') || isempty(slave.kdtree)
        slave.kdtree = KDTreeSearcher(slave.V);
    end

    % Triangle centroids on slave surface (tibia)
    F_s = slave.F;
    V_s = slave.V;
    c = (V_s(F_s(:,1),:) + V_s(F_s(:,2),:) + V_s(F_s(:,3),:)) / 3;  % nTri x 3

    % XY positions of slave centroids
    cXY = c(:,1:2);                    % nTri x 2

    % KD-tree in XY plane for master (femur)
    masterXYtree = KDTreeSearcher(master.V(:,1:2));

    % For each tibial centroid, find femur vertex with closest (x,y)
    [idxXY, ~] = knnsearch(masterXYtree, cXY);

    % Vertical gap: femur_z - tibia_z
    z_tib = c(:,3);                    % slave centroid z
    z_fem = master.V(idxXY,3);         % matched femur vertex z
    dz    = z_fem - z_tib;             % positive = femur above tibia

    % Minimum positive vertical gap, ignore overlaps
    dz_pos = dz(dz > 0);
    if isempty(dz_pos)
        e_min = 0;
    else
        e_min = min(dz_pos);
    end

    % Effective tolerance = vertical gap + baseTol
    tol_eff = e_min + baseTol;
    [A_contact, contactMask] = computeContactArea_STS(slave, master, tol_eff, opts);
end
