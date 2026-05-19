function self = set_flexion_min(self, angle)
    % Convenience function to set the minimum angle of a flexion arc
    % i.e. lowest angle in the arc is set to ANGLE, all others are relative to it.
    % Defaults to 0 degrees.
    % Equivalent to SET_ANGLE(ANGLE, @MIN);
    arguments
        self Trajectory
        angle (1,1) {mustBeNumeric} = 0
    end

    self.set_angle(angle, @min);
end
