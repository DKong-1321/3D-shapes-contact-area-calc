classdef SPSS
    properties
        Metadata
        Kinematics
    end
    methods
        function self = SPSS(envelope, interval)
            arguments
                envelope Envelope
                interval
            end

            states = envelope.states;
            directions = envelope.directions;
            specimens = envelope.specimens;
            signals = envelope.signals;


            for sg = 1:numel(signals)
                signal = signals(sg);
                for d = 1:numel(directions)
                    direction = directions(d);

                    % Preallocation stuff
                        d_prealoc = envelope.Kinematics.(specimens(1)).(states(1)).(direction).(signal);
                        ang = 1:interval:height(d_prealoc);
                        d_prealoc = d_prealoc(ang, :);
                        headers = d_prealoc.Properties.VariableNames;
                        val = nan(numel(ang), numel(states), numel(headers));
                        all_specimens = nan(numel(specimens), height(states) * numel(ang), numel(headers));
                    %
                    for sp = 1:numel(specimens)
                        specimen = specimens(sp);
                        for st = 1:numel(states)
                            state = states(st);
                            datum = envelope.Kinematics.(specimen).(state).(direction).(signal);
                            val(:, st, :) = datum{ang, :};
                        end
                        res = reshape(pagetranspose(val), [], 1, numel(headers));
                        all_specimens(sp, :, :) = pagetranspose(res);
                    end
                    for h = 1:numel(headers)
                        header = headers{h};
                        self.Kinematics.(signal).(direction).(header) = all_specimens(:, :, h);
                    end

                    % Precompute sizes
                    nA = numel(ang);
                    nS = numel(states);

                    self.Metadata.(signal).(direction).states.groups = repmat((0:nS-1)', nA, 1)';
                    self.Metadata.(signal).(direction).states.values = repmat(states(:), nA, 1)';
                    self.Metadata.(signal).(direction).angles.groups = repelem((0:nA-1)', nS)';
                    self.Metadata.(signal).(direction).angles.values = repelem((ang(:) - 1), nS)';
                end
            end
        end

        function print_to_file(self, root)
            signals = fields(self.Kinematics);
            for s = 1:numel(signals)
                signal = signals{s};
                directions = fields(self.Kinematics.(signal));
                for d = 1:numel(directions)
                    direction = directions{d};
                    headers = fields(self.Kinematics.(signal).(direction));
                    for h = 1:numel(headers)
                        header = headers{h};
                        data = self.Kinematics.(signal).(direction).(header);

                        path = fullfile(root, 'results', 'spss', signal, direction);
                        if ~exist(path, "dir")
                            mkdir(path)
                        end

                        writematrix(data, fullfile(path, [header, '.csv']));
                    end
                    metadata = self.Metadata.(signal).(direction);
                    path_metadata = fullfile(root, 'results', 'spss_metadata', signal, direction);
                    if ~exist(path_metadata, "dir")
                        mkdir(path_metadata)
                    end

                    writematrix(metadata.states.groups, fullfile(path_metadata, [header, '_states_groups.csv']));
                    writematrix(metadata.states.values, fullfile(path_metadata, [header, '_states.csv']));

                    writematrix(metadata.angles.groups, fullfile(path_metadata, [header, '_angles_groups.csv']));
                    writematrix(metadata.angles.values, fullfile(path_metadata, [header, '_angles.csv']));

                end

            end
        end
    end
end
