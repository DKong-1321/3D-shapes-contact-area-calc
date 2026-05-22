function flats = flatten(self)
    arguments
        self Path
    end
    signals = fields(self.Kinematics);
    for sg = 1:numel(signals)
        i = 1;
        signal = signals{sg};
        states = fields(self.Kinematics.(signal));
        for st = 1:numel(states)
            state = states{st};
            loading_conditions = fields(self.Kinematics.(signal).(state));
            for lc = 1:numel(loading_conditions)
                loading_condition = loading_conditions{lc};
                specimens = fields(self.Kinematics.(signal).(state).(loading_condition));
                for sp = 1:numel(specimens)
                    specimen = specimens{sp};


                    flat.state = categorical(string(state));
                    flat.loading_condition = categorical(string(loading_condition));
                    flat.specimen = categorical(string(specimen));
                    flat.kinematics = self.Kinematics.(signal).(state).(loading_condition).(specimen);

                    flats.(signal)(i) = flat;
                    i = i + 1;
                end
            end
        end
    end

end
