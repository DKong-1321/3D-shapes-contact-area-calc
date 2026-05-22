function plots = plot(self, orientations)
arguments
    self PathAverage
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
        loading_condition = loading_conditions(lc);
        figure;
        means = [];
        for s = 1:numel(states)
            state = states(s);

            if ~ismember(state, fields(self.Kinematics.(signal)))
                continue
            end

            if isempty(orientations)
                DOFs = self.Kinematics.(signal).(state).(loading_condition).mean.Properties.VariableNames;
                DOFs = setdiff(DOFs, 'flexion');
            end
            for o = 1:numel(DOFs)
                dof = DOFs{o};
                y = self.Kinematics.(signal).(state).(loading_condition).mean.(dof);
                means(s, o) = mean(y, "omitmissing");
            end
        end

        means(means == 0) = NaN;
        [~, idx_top] = max(means);
        [~, idx_bottom] = min(means);
        
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
                p = plot(x, y, 'Color', colour);


                idx = 1:10+2*s:numel(x);
                y_std = self.Kinematics.(signal).(state).(loading_condition).std.(dof);
                if s == idx_top(:, o)
                    errorbar(x(idx), y(idx), 0, y_std(idx), 'LineStyle', 'none', 'Color', colour*0.7);
                elseif s == idx_bottom(:, o)
                    errorbar(x(idx), y(idx), y_std(idx), 0,  'LineStyle', 'none', 'Color', colour*0.7);
                end

                grid on;
                axis square;
                xlabel("Flexion angle");
                ylabel(replace(dof, '_', ' '));
                title(['Loading condition: ' replace(loading_condition, '_', ' ')]);
            end
            plots.(signal).(loading_condition)(s) = p;

        end
        % sgtitle(replace(loading_condition, '_', ' '));

        for o = 1:numel(DOFs)
            nexttile(o);
            has_data = ~arrayfun(@(x) isa(x, 'matlab.graphics.GraphicsPlaceholder'), plots.(signal).(loading_condition));
            legend(plots.(signal).(loading_condition)(has_data), state_regex_inv(states(has_data)));

            path = fullfile(self.Root, 'results', 'plots', loading_condition);
            mkdir(path);
            exportgraphics(ax(o), fullfile(path, [DOFs{o} '.svg']) ,"ContentType", "vector")
        end
    end

end
end
