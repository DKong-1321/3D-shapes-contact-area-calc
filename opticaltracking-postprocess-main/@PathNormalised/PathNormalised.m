classdef PathNormalised
    properties
        Specimens
        States
        LoadingCondition
        Kinematics
        Signals
        Root
        NativeStdev
        NameNeutral
        NameNative
        NameFallback
    end
    methods
        function self = PathNormalised(path, neutral, native, fallback)
            arguments
                path
                neutral
                native
                fallback
            end

            self.Specimens = path.Specimens;
            self.States = path.States;
            self.LoadingCondition = path.LoadingCondition;
            self.Kinematics = path.Kinematics;
            self.Signals = path.Signals;
            self.Root = path.Root;
            self.NameNeutral = neutral;
            self.NameNative = native;
            self.NameFallback = fallback;

            signals = [path.Signals];


            for sg = 1:numel(signals)
                signal = signals(sg);

                if ~ismember(signal, fields(path.Kinematics))
                    continue
                end
                datum = path.Kinematics.(signal);

                states = string(fields(datum));
                if ismember(native, states)
                    native_neutral = datum.(native).(neutral);
                else
                    if ismember(fallback, states)
                        native_neutral = datum.(fallback).(neutral);
                    else
                        error("Couldn't find native neutral or use fallback");
                    end
                end


                for st = 1:numel(states)
                    state = states(st);
                    loading_conditions = string(fields(datum.(state)));

                    for lc = 1:numel(loading_conditions)
                        loading_condition = loading_conditions(lc);
                        specimens = string(fields(datum.(state).(loading_condition)));


                        for sp = 1:numel(specimens)
                            specimen = specimens(sp);

                            current = datum.(state).(loading_condition).(specimen);
                            if isempty(current)
                                continue
                            end

                            current_native_neutral = native_neutral.(specimen);

                            quantised = quantise({current, current_native_neutral});

                            quantised_current = fillmissing(quantised{1}, "makima", "EndValue", "none");
                            quantised_native_neutral = fillmissing(quantised{2}, "makima", "EndValue", "none");

                            self.Kinematics.(signal).(state).(loading_condition).(specimen) = quantised_current - quantised_native_neutral;
                            self.Kinematics.(signal).(state).(loading_condition).(specimen).flexion = quantised_current.flexion;
                        end
                    end

                end
                for sp = 1:numel(specimens)
                    specimen = specimens(sp);
                    a{sp} = native_neutral.(specimen);
                end
                if all(cellfun(@isempty, a))
                    continue
                end
                tables = cellfun(@table2array, quantise(a), 'UniformOutput', false);
                stacked = cat(3, tables{:});
                stacked = fillmissing(stacked, "pchip", "EndValues", "none");
                stdev = std(stacked, 0, 3);

                self.NativeStdev.(signal) = stdev;
            end
        end

        function [neutral, native, fallback] = names(self)
            neutral = self.NameNeutral;
            native = self.NameNative;
            fallback = self.NameFallback;
        end
    end

end
