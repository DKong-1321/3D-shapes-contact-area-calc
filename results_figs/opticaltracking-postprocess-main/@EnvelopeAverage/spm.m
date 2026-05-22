function spmi = spm(self)
    states = self.States;
    directions = self.Directions;
    signals = self.Signals;

    for sg = 1:numel(signals)
        signal = signals(sg);
        for d = 1:numel(directions)
            direction = directions(d);
            x = self.Kinematics.(states(1)).(signal).(direction).mean;
            val = nan(height(x), numel(states), width(x));

            for st = 1:numel(states)
                state = states(st);
                datum = self.Kinematics.(state).(signal).(direction).mean;
                val(:, st, :) = table2array(datum);
            end

            headers = datum.Properties.VariableNames;
            for h = 1:numel(headers)
                header = headers{h};
                spm = spm1d.stats.anova1rm(val(:, :, h)', 1:numel(states));
                spmi.(signal).(direction).(header) = spm.inference(0.05);
            end
        end
    end
end
