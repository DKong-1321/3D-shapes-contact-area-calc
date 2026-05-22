function self = intraspecimen_mean(self)
    for t = 1:numel(self)
        data = self(t);
        signals = data.signals;
        for sg = 1:numel(signals)
            signal = signals(sg);
            datum = data.Kinematics.(signal);
            if isempty(datum)
                continue
            end

            self(t).Kinematics.(signal) = calc_means(self(t).Kinematics.(signal));
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
            self(t).Sensors.(sensor_name) = calc_means(sensor);
        end

    end
end

function means = calc_means(datum)
    headers = datum.Properties.VariableNames;
    is_flexion = strcmpi(headers, 'flexion');
    flexion = headers{is_flexion};

    minima = find_minima(datum.(flexion));
    runs = split_run(datum, minima);
    if isempty(runs)
        means = [];
        return
    end
    quantised_runs = quantise(runs);
    mat = cellfun(@table2array, quantised_runs, "UniformOutput", false);
    mat_stack = cat(3, mat{:});

    mat_means = mean(mat_stack, 3, "omitmissing");
    tab_means = array2table(mat_means, "VariableNames",headers);

    means = fillmissing(tab_means, "pchip");
end
