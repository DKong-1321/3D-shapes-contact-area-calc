function q = angle2quaternion(roll, pitch, yaw)
% ZYX order. Implementation copied from https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
cr = cosd(roll/2);
sr = sind(roll/2);
cp = cosd(pitch/2);
sp = sind(pitch/2);
cy = cosd(yaw/2);
sy = sind(yaw/2);

qw = cr * cp * cy + sr * sp * sy;
qx = sr * cp * cy - cr * sp * sy;
qy = cr * sp * cy + sr * cp * sy;
qz = cr * cp * sy - sr * sp * cy;

q = [qw qx qy qz];

end