function plots = plot(self, orientations)
arguments
    self PathNormalisedAverage
    orientations = []
end
signals = self.Signals;
states = self.States;
colours = lines(numel(states));
loading_conditions = [self.LoadingCondition];

for sg = 1:numel(signals)
    signal = signals(sg);

    if ~ismember(signal, fields(self.Kinematics))
        continue
    end

    for lc = 1:numel(loading_conditions)
        figure;
        loading_condition = loading_conditions(lc);

        [idx_top, idx_bottom] = find_top_bottom(self.Kinematics.(signal),  loading_condition, states, orientations, self.NameNative);

        for s = 1:numel(states)
            state = states(s);


            colour = colours(s, :);

            if ~ismember(state, fields(self.Kinematics.(signal)))
                continue
            end

            if isempty(orientations)
                DOFs = self.Kinematics.(signal).(state).(loading_condition).mean.Properties.VariableNames;
                DOFs = setdiff(DOFs, 'flexion');
            end

            for o = 1:numel(DOFs)
                dof = DOFs{o};
                ax(o) = nexttile(o); hold on;
                x = self.Kinematics.(signal).(state).(loading_condition).mean.flexion;
                y = self.Kinematics.(signal).(state).(loading_condition).mean.(dof);
                y = smoothdata(y, "gaussian", 5);

                if state == self.NameNative
                    y_std = self.NativeStdev.(signal);
                    len = min(size(y, 1), size(y_std, 1));
                    y_std = y_std(1:len, o);
                    y = y(1:len);
                    x = x(1:len);
                    y_upper = y + y_std;
                    y_lower = y - y_std;

                    peak = find(x == max(x), 1, 'first');
                    segments = {1:peak, peak:numel(x)};

                    for k = 1:numel(segments)
                        idx = segments{k};
                        fill([x(idx); flipud(x(idx))], [y_upper(idx); flipud(y_lower(idx))], colour, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    end
                    if loading_condition ~= self.NameNeutral
                        p = plot(x, y, 'Color', colour);
                    end
                    % NameFallback 
                else
                    p = plot(x, y, 'Color', colour);

                    idx = 1:10+2*s:numel(x);
                    y_std = self.Kinematics.(signal).(state).(loading_condition).std.(dof);
                    if s == idx_top(:, o)
                        errorbar(x(idx), y(idx), 0, y_std(idx), 'LineStyle', 'none', 'Color', colour*0.7);
                    elseif s == idx_bottom(:, o)
                        errorbar(x(idx), y(idx), y_std(idx), 0,  'LineStyle', 'none', 'Color', colour*0.7);
                    end
                end



                grid on;
                axis square;
                xlabel("Flexion angle");
                ylabel(replace(dof, '_', ' '));
            end
            plots.(signal).(loading_condition)(s) = p;

        end
        % sgtitle(replace(loading_condition, '_', ' '));

        for o = 1:numel(DOFs)
            nexttile(o);
            has_data = ~arrayfun(@(x) isa(x, 'matlab.graphics.GraphicsPlaceholder'), plots.(signal).(loading_condition));
            legend(plots.(signal).(loading_condition)(has_data), state_regex_inv(states(has_data)));
            title(['Normalised path' ['Loading condition: ', replace(loading_condition, '_', ' ')]]);

            path = fullfile(self.Root, 'results', 'plots', loading_condition);
            mkdir(path);
            exportgraphics(ax(o), fullfile(path, [DOFs{o} '.svg']) ,"ContentType", "vector")
        end
    end

end
end

function [idx_top, idx_bottom] = find_top_bottom(data, loading_condition, states, orientations, native)

    means = [];
    for s = 1:numel(states)
        state = states(s);

        if ~ismember(state, fields(data))
            continue
        end


        if isempty(orientations)
            DOFs = data.(state).(loading_condition).mean.Properties.VariableNames;
            DOFs = setdiff(DOFs, 'flexion');
        end
        for o = 1:numel(DOFs)
            dof = DOFs{o};
            y = data.(state).(loading_condition).mean.(dof);
            means(s, o) = mean(y, "omitmissing");
        end
    end

    means(means == 0) = NaN;
    means(states == native, :) = NaN;
    [~, idx_top] = max(means);
    [~, idx_bottom] = min(means);
end
