function problems = check_health(self, min_flexion_arc, fraction_of_data)
% Checks all data for a minimum flexion arc and missing data
% Defaults to a minimum flexion arc of 60 degrees and 20% of data missing.
    arguments
        self Trajectory
        min_flexion_arc = 60
        fraction_of_data = 0.2;
    end
    i = 1;
    for t = 1:numel(self)
        data = self(t);
        signals = data.signals;

        for sg = 1:numel(signals)
            signal = signals(sg);
            datum = data.Kinematics.(signal);

            if isempty(datum)
                continue
            end

            headers = datum.Properties.VariableNames;
            is_flexion = strcmpi(headers, 'flexion');
            flexion = headers{is_flexion};

            arc = abs(max(datum.(flexion)) - min(datum.(flexion)));
            has_min_flex = arc > min_flexion_arc;

            num_nan = sum(isnan(table2array(datum)), "all");
            complete_data = num_nan < fraction_of_data * numel(datum);

            if has_min_flex && complete_data
                continue
            end

            problems(i).SpecimenName = data.SpecimenName;
            problems(i).SpecimenState = data.SpecimenState;
            problems(i).LoadingCondition = data.LoadingCondition;
            problems(i).(signal).flexion_arc = arc;
            problems(i).(signal).missing_data = num_nan/height(datum);
            i = i + 1;
        end
    end
end
