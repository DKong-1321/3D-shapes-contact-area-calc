function self = print_to_file(self, path)
    % % Prepare the folders
    fp_results = fullfile(path, "results", "average", "stability_envelope");
    states = self.States;
    signals = self.Signals;
    directions = self.Directions;
    for st = 1:numel(states)
        state = states(st);
        for sg = 1:numel(signals)
            signal = signals(sg);

            filepath = fullfile(fp_results, signal, state);
            mkdir(filepath);
            for d = 1:numel(directions)
                direction = directions(d);

                datum = self.Kinematics.(state).(signal).(direction);
                headers = datum.mean.Properties.VariableNames;

                t_mean = datum.mean;
                t_std = datum.std;

                t_mean.Properties.VariableNames = headers + "_mean";
                t_std.Properties.VariableNames = headers + "_std";


                writetable([t_mean t_std], strcat(fullfile(filepath, direction), '.csv'));
            end
        end
    end
end
