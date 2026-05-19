function self = print_to_file(self, path)
% % Prepare the folders
fp_results = fullfile(path, "results");
states = unique(self.states);
signals = unique(self.signals);
specimens = unique(self.specimen);
loading_conditions = unique([self.LoadingCondition]);
for sg = 1:numel(signals)
    signal = signals(sg);
    for sp = 1:numel(specimens)
        specimen = specimens(sp);
        for st = 1:numel(states)
            state = states(st);
            for lc = 1:numel(loading_conditions)
                loading_condition = loading_conditions(lc);

                is_specimen = specimen == [self.specimen];
                is_state = state == [self.states];
                is_lc = loading_condition == [self.LoadingCondition];


                is_current = is_specimen & is_state & is_lc;
                current_specimen = [self(is_current).Kinematics];
                for n = 1:numel(current_specimen)
                    if n > 1
                        keyboard
                    end
                    if ~ismember(signal, fields(current_specimen))
                        continue
                    end
                    datum = current_specimen(n).(signal);
                    if isempty(datum)
                        continue
                    end

                    filepath = fullfile(fp_results, signal, state, loading_condition);
                    mkdir(filepath);

                    writetable(datum, strcat(fullfile(filepath, specimen), '.csv'));
                end
            end
        end
    end
end
end
