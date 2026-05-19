function [Rx, Ry, Rz] = quaternion2euler(qw,qx,qy,qz)
    %Converts quaternions to euler angles: http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    Rx=atan2((2*(qw .* qx + qy .* qz)), (1 - 2*(qx.^2+qy.^2)))*180/pi;
    Ry=asind(2*(qw.*qy-qz.*qx));
    Rz=atan2((2*(qw.*qz+qx.*qy)),(1-2*(qy.^2+qz.^2)))*180/pi;
end

