function [flex, ext] = split_flex_ext(self)
    states = self.States;
    signals = self.Signals;
    specimens = self.Specimens;
    loading_conditions = self.LoadingCondition;

    flex = self;
    ext = self;

    flex.Kinematics = [];
    ext.Kinematics = [];

    for sg = 1:numel(signals)
        signal = signals(sg);
        for st = 1:numel(states)
            state = states(st);
            for lc = 1:numel(loading_conditions)
                loading_condition = loading_conditions(lc);
                for sp = 1:numel(specimens)
                    specimen = specimens(sp);

                    curr_specimens = fields(self.Kinematics.(signal).(state).(loading_condition));
                    if ~ismember(specimen, curr_specimens)
                        continue
                    end
                    datum = self.Kinematics.(signal).(state).(loading_condition).(specimen);
                    if isempty(datum)
                        continue
                    end

                    [~, n] = max(datum.flexion);

                    flex.Kinematics.(signal).(state).(loading_condition).(specimen) = datum(1:n, :);
                    ext.Kinematics.(signal).(state).(loading_condition).(specimen) = datum(n:end, :);
                end
            end
        end
    end
end
