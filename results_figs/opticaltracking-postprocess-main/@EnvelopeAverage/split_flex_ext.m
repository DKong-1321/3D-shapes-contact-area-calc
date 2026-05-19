function [flex, ext] = split_flex_ext(self)
directions = self.Directions;
states = self.States;
signals = self.Signals;

flex = self;
ext = self;

for d = 1:numel(directions)
    direction = directions(d);
    for st = 1:numel(states)
        state = states(st);

        for sg = 1:numel(signals)
            signal = signals(sg);
            datum = self.Kinematics.(state).(signal).(direction);
            headers = fieldnames(datum);

            is_flexion = strcmpi(headers, 'flexion');
            flexion = headers{is_flexion};
            [~, n] = max(flexion);
            flex.Kinematics.(state).(signal).(direction) = datum.(header)(1:n, :);
            ext.Kinematics.(state).(signal).(direction) = datum.(header)(n:end, :);
        end
    end
end
end
