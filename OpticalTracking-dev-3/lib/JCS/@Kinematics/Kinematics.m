classdef Kinematics
    properties
        trajectories
    end
    properties(Access = private)
        digitisation
        config
    end
    methods
        function self = Kinematics(trajectories, digitisation, config)
            self.trajectories = trajectories;
            self.config = config;
            self.digitisation = digitisation;
        end
    end
end

