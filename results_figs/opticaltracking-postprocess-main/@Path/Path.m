classdef Path
    properties
        Specimens
        States
        LoadingCondition
        Kinematics
        Signals
        Root
    end

    methods
        function self = Path(data, specimens, states, loading_conditions)
            self.Specimens = specimens;
            self.States = setdiff(unique(states), ["UKA_w_pACL", "Unoptimised"]);
            self.LoadingCondition = loading_conditions;
            self.Signals = data.signals;
            self.Root = data.Root;

            for sp = 1:numel(specimens)
                specimen = specimens(sp);
                is_specimen = [data.SpecimenName] == specimen;
                for st = 1:numel(states)
                    state = states(st);
                    is_state = [data.SpecimenState] == state;
                    for lc = 1:numel(loading_conditions)

                        loading_condition = loading_conditions(lc);
                        is_loading_condition = [data.LoadingCondition] == loading_condition;

                        mask = is_specimen & is_loading_condition & is_state;
                        datum = [data(mask)];
                        datum = datum(~isempty(datum));
                        if isempty(datum)
                            continue
                        end
                        signals = fields([datum.Kinematics]);
                        for sg = 1:numel(signals)
                            signal = signals{sg};
                            try
                            self.Kinematics.(signal).(state).(loading_condition).(specimen) = datum.Kinematics.(signal);
                            catch me
                                keyboard
                            end
                        end
                    end
                end
            end
        end
        function data = neutral(self)
            data = self;
            states = data.States;
            for st = 1:numel(states)
                state = states(st);

                datum = [data.Kinematics.(state)];
                is_neutral = contains(fields(datum), "neutral", "IgnoreCase", true);
                not_neutrals = fields(datum);
                not_neutrals = not_neutrals(~is_neutral);
                for n = 1:numel(not_neutrals)
                    not_neutral = not_neutrals{n};
                    data.Kinematics.(state).(not_neutral) = [];
                end
            end
        end
        function o = average(self)
            o = PathAverage(self);
        end
        function o = filter_signal(self, signal)
            self.Signals = self.Signals(contains(self.Signals, signal));
            o = self;
        end
        function o = exclude_specimen(self, specimen)
            o = self;
            mask = contains(o.Specimens, specimen, "IgnoreCase", true);
            o.Specimens = o.Specimens(~mask);
        end

        function o = exclude_specimen_exact(self, specimen)
            o = self;
            specimens_remaining = setdiff(o.Specimens, specimen);
            o.Specimens = specimens_remaining;
        end

    end
end

