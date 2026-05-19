classdef PathAverage
    properties
        Specimens
        States
        LoadingCondition
        Signals
        Kinematics
        Root
    end
    methods
        function self = PathAverage(paths)
            arguments
                paths
            end

            states = paths.States;
            signals = paths.Signals;
            specimens = paths.Specimens;
            loading_conditions = paths.LoadingCondition;

            for sg = 1:numel(signals)
                signal = signals(sg);
                if ~ismember(signal, fields(paths.Kinematics))
                    continue
                end
                for st = 1:numel(states)
                    state = states(st);
                    if ~ismember(state, fields(paths.Kinematics.(signal)))
                        continue
                    end

                    for lc = 1:numel(loading_conditions)
                        loading_condition = loading_conditions(lc);
                        data = [paths.Kinematics.(signal).(state).(loading_condition)];
                        tables = cell(size(specimens));
                        for sp = 1:numel(specimens)
                            specimen = specimens(sp);
                            tables{sp} = data.(specimen);
                        end

                        is_empty = cellfun(@isempty, tables);
                        tables = tables(~is_empty);

                        if isempty(tables)
                            continue
                        end

                        tables = quantise(tables);

                        tables = cellfun(@table2array, tables, "UniformOutput", false);
                        stacked = cat(3, tables{:});
                        % if any(isnan(stacked), "all")
                        %     keyboard
                        % end
                        stacked = fillmissing(stacked, "pchip", "EndValues", "none");
                        avg = mean(stacked, 3);
                        stdev = std(stacked, 0, 3);

                        headers = data.(specimen).Properties.VariableNames;

                        self.Kinematics.(signal).(state).(loading_condition).mean = array2table(avg, "VariableNames", headers);
                        self.Kinematics.(signal).(state).(loading_condition).std = array2table(stdev, "VariableNames", headers);
                    end
                end
            end
            self.LoadingCondition = paths.LoadingCondition;
            self.Signals = paths.Signals;
            self.States = paths.States;
            self.Root = paths.Root;
        end
    end
end
