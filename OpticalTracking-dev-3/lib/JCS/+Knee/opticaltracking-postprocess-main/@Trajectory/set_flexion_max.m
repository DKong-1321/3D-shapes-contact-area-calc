function self = set_flexion_max(self, angle)
    % Convenience function to set the maximum angle of a flexion arc
    % Defaults to 90 degrees.
    % Equivalent to SET_ANGLE(ANGLE, @MAX);
    arguments
        self Trajectory
        angle (1,1) {mustBeNumeric} = 90
    end

    self.set_angle(angle, @max);
end
