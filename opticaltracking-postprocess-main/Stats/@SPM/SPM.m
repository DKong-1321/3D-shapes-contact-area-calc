classdef SPM
    properties
        % rawdata
        % stats
        signal
        dof
        within_subject
        between_subject
        states
        specimens
        loading_conditions
    end

    methods
        function self = SPM(data)
            if nargin > 0
                signals = fields(data);
                total_iter = 0;
                for sg = 1:numel(signals)
                    signal = signals{sg};

                    if isempty(data.(signal))
                        continue
                    end
                    datum = data.(signal);
                    headers = datum(1).kinematics.Properties.VariableNames;
                    total_iter = total_iter + numel(headers);
                end

                self(total_iter).signal = [];
                self(total_iter).dof = [];
                self(total_iter).between_subject = [];
                self(total_iter).within_subject = [];

                k = 1;
                for sg = 1:numel(signals)
                    signal = signals{sg};

                    if isempty(data.(signal))
                        continue
                    end

                    states = findgroups([data.(signal).state]);
                    loading_conditions = findgroups([data.(signal).loading_condition]);
                    specimens = findgroups([data.(signal).specimen]);

                    datum = data.(signal);


                    headers = datum(1).kinematics.Properties.VariableNames;

                    quantised = quantise({datum.kinematics});

                    h = min(cellfun(@height, quantised));
                    val = nan(h, numel(unique(specimens)) * numel(unique(states)) * numel(unique(loading_conditions)), numel(headers));

                    for i = 1:numel(datum)
                        arr = table2array(quantised{i});
                        val(1:size(arr,1), i, 1:size(arr, 2)) = arr;
                    end

                    val = fillmissing(val, "makima", "EndValues", "none");

                    % self.rawdata.(signal) = datum;
                    for h = 1:numel(headers)
                        header = headers{h};
                        self(k).signal = categorical({signal});
                        self(k).dof = categorical({header});
                        self(k).states = categories([datum.state]);
                        self(k).specimens = categories([datum.specimen]);
                        self(k).loading_conditions = categories([datum.loading_condition]);
                        self(k).between_subject = spm1d.stats.anova2(val(:, :, h)', states, loading_conditions);
                        self(k).within_subject = spm1d.stats.anova2rm(val(:, :, h)', states, loading_conditions, specimens);
                        
                        % self.between_subject.(signal).(header) = spm1d.stats.anova2(val(:, :, h)', states, loading_conditions);
                        % self.within_subject.(signal).(header) = spm1d.stats.anova2rm(val(:, :, h)', states, loading_conditions, specimens);

                        k = k + 1;
                    end

                end
            end
        end

        function o = dunnett(self, control)
            if nargin < 2
                error("Missing control group. options: %s", strjoin(self.states, ', '))
            end
            if ~ismember(control, self.states)
                error("Not a valid control. options: %s", strjoin(self.states, ', '))
            end

            alpha = 0.05;
            n_tests = numel(self.states) - 1;

            p_critical = spm1d.util.p_corrected_bonf(alpha, n_tests);


            o = Dunnett(self, control, p_critical);
        end

    end
end
