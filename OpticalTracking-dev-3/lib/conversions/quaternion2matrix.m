function R = quaternion2matrix(q)
% Based on https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

qw = q(:, 1);
qx = q(:, 2);
qy = q(:, 3);
qz = q(:, 4);
if norm(q) ~= 1 % Normalise
    q = q/norm(q);
end

% % Inhomogeneous
% R = [
%       1-2*(qy^2 + qz^2),   2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw),    0;
%       2*(qx*qy + qz*qw),   1-2*(qx^2 + qz^2),     2*(qy*qz - qx*qw),    0;
%       2*(qx*qz - qy*qw),   2*(qy*qz + qx*qw),     1 - 2*(qx^2 + qy^2),  0;
%       0,                       0,                        0,             1;
%     ];

% Homogeneous
R = repmat(eye(4), [1, 1, numel(qw)]);
R(1,1,:) = qw.^2 + qx.^2 - qy.^2 -qz.^2;
R(2,1,:) = 2*(qx.*qy + qw.*qz);
R(3,1,:) = 2*(qx.*qz - qw.*qy);

R(1,2,:) = 2*(qx.*qy - qw.*qz);
R(2,2,:) = qw.^2 - qx.^2 + qy.^2 - qz.^2;
R(3,2,:) = 2*(qw.*qx + qy.*qz);

R(1,3,:) = 2*(qw.*qy + qx.*qz);
R(2,3,:) = 2*(qy.*qz - qw.*qx);
R(3,3,:) = qw.^2 - qx.^2 - qy.^2 + qz.^2;

end