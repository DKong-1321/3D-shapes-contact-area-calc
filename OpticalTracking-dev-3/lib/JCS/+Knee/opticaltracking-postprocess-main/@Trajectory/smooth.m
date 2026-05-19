function self = smooth(self, smoothing_func, fillmissing_func)
% Fills missing data and applies smoothing function to all trajectories. If no smoothing function is provided,
% the default @(x) smoothdata(x, "gaussian", 10); is used.
% If no fillmissing_func is provided, uses @(x) fillmissing(x, 'pchip")
% SMOOTHING_FUNC, FILLMISSING_FUNC must be an anonymous function.
arguments
    self Trajectory
    smoothing_func = @(x) smoothdata(x, "gaussian", 10);
    fillmissing_func = @(x) fillmissing(x, "pchip"); % Missing quanta are filled with this function. E.g., have flexion angles 40 and 42; it fills the 41 linearly.
end
for t = 1:numel(self)
    data = self(t);
    signals = data.signals;
    for sg = 1:numel(signals)
        signal = signals(sg);
        datum = data.Kinematics.(signal);
        if isempty(datum)
            continue
        end
        headers = datum.Properties.VariableNames;
        R = table2array(datum);
        R(any(R==0, 2),:) = NaN;
        R = fillmissing_func(R);
        R = smoothing_func(R);

        self(t).Kinematics.(signal) = array2table(R, "VariableNames", headers);
    end
end
end
