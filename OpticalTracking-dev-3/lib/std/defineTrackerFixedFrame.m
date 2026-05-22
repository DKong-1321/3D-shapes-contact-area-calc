function [ gTt0 ] = defineTrackerFixedFrame(eulerAngles,XYZ)
    if isempty(eulerAngles) || isempty(XYZ)
        gTt0 = [];
        return;
    end 
    identity_mat = @(n) repmat(eye(4), 1, 1, n);
    n = height(eulerAngles);

    Qx = identity_mat(n);
    Qy = identity_mat(n);
    Qz = identity_mat(n);
    trans = identity_mat(n);

    Rx = eulerAngles(:, 1);
    Ry = eulerAngles(:, 2);
    Rz = eulerAngles(:, 3);

    Qx(2:3, 2:3, :) = [cosd(Rx), -sind(Rx); sind(Rx), cosd(Rx)];

    Qy(1, 1, :) = cosd(Ry);
    Qy(1, 3, :) = sind(Ry);
    Qy(3, 1, :) = -sind(Ry);
    Qy(3, 3, :) = cosd(Ry);

    Qz(1:2, 1:2, :) = [cosd(Rz), -sind(Rz); sind(Rz), cosd(Rz)];

    rot = pagemtimes(pagemtimes(Qz, Qy), Qx);

    trans(1:3, 4, :) = XYZ';

    gTt0 = pagemtimes(trans, rot);
end
