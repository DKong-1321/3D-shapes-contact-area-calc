function [minimum, maximum] = flexion_extrema(self)
    maximum = -inf;
    minimum = inf;

    data = [self.Kinematics];
    signals = self.signals;
    for sg = 1:numel(signals)
        signal = signals(sg);
        datum = {data.(signal)};
        is_empty = cellfun(@isempty, datum);
        datum = datum(~is_empty);

        if isempty(datum)
            continue
        end

        local_maxima = cellfun(@(x) max(x.flexion), datum, "UniformOutput", false);
        local_maximum = max([local_maxima{:}]);
        local_minima = cellfun(@(x) min(x.flexion), datum, "UniformOutput", false);
        local_minimum = min([local_minima{:}]);

        maximum = max(local_maximum, maximum);
        minimum = min(local_minimum, minimum);
    end
end
