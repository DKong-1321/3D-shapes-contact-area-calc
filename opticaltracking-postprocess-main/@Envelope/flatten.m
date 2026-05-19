function flats = flatten(self)
    arguments
        self Envelope
    end
    specimens = fields(self.Kinematics);
    for sp = 1:numel(specimens)
        i = 1;
        specimen = specimens{sp};
        states = fields(self.Kinematics.(specimen));
        for st = 1:numel(states)
            state = states{st};
            loading_conditions = fields(self.Kinematics.(specimen).(state));
            for lc = 1:numel(loading_conditions)
                loading_condition = loading_conditions{lc};
                
                datum = self.Kinematics.(specimen).(state).(loading_condition);
                if isempty(datum)
                    continue
                end
                signals = fields(datum);
                for sg = 1:numel(signals)
                    signal = signals{sg};


                    flat.state = categorical(string(state));
                    flat.loading_condition = categorical(string(loading_condition));
                    flat.specimen = categorical(string(specimen));
                    flat.kinematics = datum.(signal);

                    flats.(signal)(i) = flat;
                    i = i + 1;
                end
            end
        end
    end

    % for sg = 1:numel(signals)
    %     signal = signals{sg};
    % is_empty = cellfun(@isempty, {flats.(signal).kinematics});
    % out.(signal) = flats.(signal)(~is_empty);
    % end
end
