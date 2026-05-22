function self = print_to_file(self, path)
    % % Prepare the folders
    fp_results = fullfile(path, "results", "per_specimen", "stability_envelope");
    states = self.states;
    signals = self.signals;
    directions = self.directions;
    specimens = self.specimens;
    for sp = 1:numel(specimens)
        specimen = specimens(sp);
        for st = 1:numel(states)
            state = states(st);
            for sg = 1:numel(signals)
                signal = signals(sg);

                filepath = fullfile(fp_results, signal, state, specimen);
                mkdir(filepath);
                for d = 1:numel(directions)
                    direction = directions(d);

                    datum = self.Kinematics.(specimen).(state).(direction).(signal);

                    writetable(datum, strcat(fullfile(filepath, direction), '.csv'));
                end
            end
        end
    end
end
