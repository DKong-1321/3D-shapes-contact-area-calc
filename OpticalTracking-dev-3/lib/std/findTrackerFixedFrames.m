function [ gTtiF ] = findTrackerFixedFrames(rigid_body)
if rigid_body.is_none
    gTtiF = [];
    return
end
rigid_body = rigid_body.unwrap();
eulerAngles = rigid_body.rotations;
XYZ = rigid_body.translations;

% Creates 4 x 4 x m matrices, where each page is gTtiF (tracker in global
% frame of reference)
Rx = reshape(eulerAngles(:,1), [1, 1, length(eulerAngles(:,1))]);
Ry = reshape(eulerAngles(:,2), [1, 1, length(eulerAngles(:,2))]);
Rz = reshape(eulerAngles(:,3), [1, 1, length(eulerAngles(:,3))]);
Qx = repmat(eye(4), 1, 1, length(eulerAngles(:,1)));
Qy = repmat(eye(4), 1, 1, length(eulerAngles(:,1)));
Qz = repmat(eye(4), 1, 1, length(eulerAngles(:,1)));

% Qx=[1   0         0
%     0   cosd(Rx)  -sind(Rx)
%     0   sind(Rx)  cosd(Rx)];
Qx(2, 2, :) = cosd(Rx);
Qx(2, 3, :) = -sind(Rx);
Qx(3, 2, :) = sind(Rx);
Qx(3, 3, :) = cosd(Rx);

% Qy=[cosd(Ry)    0   sind(Ry)
%     0           1   0
%     -sind(Ry)   0   cosd(Ry)];
Qy(1, 1, :) = cosd(Ry);
Qy(1, 3, :) = sind(Ry);
Qy(3, 1, :) = -sind(Ry);
Qy(3, 3, :) = cosd(Ry);

% Qz=[cosd(Rz)    -sind(Rz)   0
%     sind(Rz)    cosd(Rz)    0
%     0           0           1];  
Qz(1, 1, :) = cosd(Rz);
Qz(1, 2, :) = -sind(Rz);
Qz(2, 1, :) = sind(Rz);
Qz(2, 2, :) = cosd(Rz);


rot = pagemtimes(pagemtimes(Qz, Qy), Qx);

trans = repmat(eye(4), 1, 1, length(XYZ(:,1)));
trans(1, 4, :) = XYZ(:,1);
trans(2, 4, :) = XYZ(:,2);
trans(3, 4, :) = XYZ(:,3);

gTtiF = pagemtimes(trans, rot);

end
