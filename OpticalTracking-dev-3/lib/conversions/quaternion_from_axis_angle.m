function [q, qc] = quaternion_from_axis_angle(axis, angle)
% https://www.youtube.com/watch?v=zjMuIxRvygQ
axis = axis(:)';
q = [cos(angle/2),  sin(angle/2) * axis];
q = q/norm(q);
qc = [q(1) -q(2:4)];
end