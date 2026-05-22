function self = set_angle(self, angle, func)
    % FUNC must be an anonymous function. E.g. if you want maximum angle to be 120 degrees,
    % the arguments are 120, @max. If you want first angle to be 0, then arguments are: 0, @(x) x(1)
    arguments
        self Trajectory
        angle (1,1) {mustBeNumeric}
        func
    end

    for t = 1:numel(self)
        data = self(t);
        signals = data.signals;
        %% Kinematics
        for sg = 1:numel(signals)
            signal = signals(sg);
            kinematics = data.Kinematics.(signal);
            if isempty(kinematics)
                continue
            end
            [offset, flexion ] = to_zero(kinematics, angle, func);
            self(t).Kinematics.(signal).(flexion) = offset;

        end

        %% Sensors
        sensors = [self(t).Sensors];
        if isempty(sensors)
            continue
        end
        sensor_names = fields(sensors);
        for s = 1:numel(sensor_names)
            sensor_name = sensor_names{s};
            sensor = sensors.(sensor_name);
            [offset, flexion] = to_zero(sensor, angle, func);
            self(t).Sensors.(sensor_name).(flexion) = offset;
        end
    end
end

function [data_offset, flexion] = to_zero(data, angle, func)
    headers = data.Properties.VariableNames;
    is_flexion = strcmpi(headers, 'flexion');
    flexion = headers{is_flexion};

    delta = angle - func(data.(flexion));
    data_offset = data.(flexion) + delta;
end
