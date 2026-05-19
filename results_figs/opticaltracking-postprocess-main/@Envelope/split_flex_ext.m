function [flex, ext] = split_flex_ext(self)
    directions = self.directions;
    states = self.states;
    specimens = self.specimens;
    self.SpecimenName = specimens;
    flex = self;
    ext = self;

    for d = 1:numel(directions)
        direction = directions(d);
        for st = 1:numel(states)
            state = states(st);

            for sp = 1:numel(specimens)
                specimen = specimens(sp);
                signals = fieldnames(self.Kinematics.(specimen).(state).(direction));
                for sg = 1:numel(signals)
                    signal = signals{sg};
                    datum = self.Kinematics.(specimen).(state).(direction).(signal);
                    headers = datum.Properties.VariableNames;
                    is_flexion = strcmpi(headers, 'flexion');
                    flexion = headers{is_flexion};
                    [~, n] = max(flexion);
                    flex.Kinematics.(specimen).(state).(direction).(signal) = datum(1:n, :);
                    ext.Kinematics.(specimen).(state).(direction).(signal) = datum(n:end, :);
                end
            end
        end
    end

end
