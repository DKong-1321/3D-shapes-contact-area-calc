function fig = plot(self, highlight_idx)
    arguments
        self
        highlight_idx = []
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
            figure;
            fig(t) = tiledlayout(round(numel(headers)/2), 2);
            sgtitle([[self(t).SpecimenName ' ' replace(self(t).SpecimenState, '_', ' ') ' ' self(t).loading_condition] signals(sg)]);
            for h = 1:numel(headers)
                % Plot data
                nexttile(h);
                header = headers{h};
                p = datum.(header);
                plot(p);
                ylabel(replace(headers{h}, '_', ' '));
                % Plot the interpolation
                if isempty(highlight_idx)
                    continue;
                end
                interpolated_points = highlight_idx(t).(signal);

                if any(interpolated_points)
                    hold on;
                    scatter(find(interpolated_points), p(interpolated_points), 10, "filled", "red")
                    hold off;
                end
            end
        end
    end
end
