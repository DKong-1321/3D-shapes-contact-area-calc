function plots = plot(self)
arguments
    self Path
end
signals = self.Signals;
states = self.States;
colours = lines(numel(states));
specimens = self.Specimens;
loading_conditions = self.LoadingCondition;

for sg = 1:numel(signals)
    signal = signals(sg);
    if ~ismember(signal, fields(self.Kinematics))
        continue
    end
    for lc = 1:numel(loading_conditions)
        loading_condition = loading_conditions(lc);

        for sp = 1:numel(specimens)
            specimen = specimens(sp);
            fig = [];
            for s = 1:numel(states)
                state = states(s);
                colour = colours(s, :);

                if ~ismember(state, fields(self.Kinematics.(signal)))
                    continue
                end
                datum = self.Kinematics.(signal).(state).(loading_condition).(specimen);
                if isempty(datum)
                    continue
                end

                if isempty(fig)
                    fig = figure;
                end

                orientations = datum.Properties.VariableNames;
                orientations = setdiff(orientations, 'flexion');
                for o = 1:numel(orientations)


                    nexttile(o); hold on;
                    x = datum.flexion;
                    y = datum.(orientations{o});
                    plots.(signal).(loading_condition).(specimen)(s) = plot(x, y, 'Color', colour);

                    grid on;
                    axis square;
                    xlabel("Flexion angle");
                    ylabel(replace(orientations{o}, '_', ' '));
                end
            end
            if ~isempty(fig)
                sgtitle([specimen loading_condition replace(signal, '_', ' ')]);
                has_data = ~arrayfun(@(o) isa(o, 'matlab.graphics.GraphicsPlaceholder'), plots.(signal).(loading_condition).(specimen));
                legend(plots.(signal).(loading_condition).(specimen)(has_data), state_regex_inv(states(has_data)));
            end
        end
    end
end
end


function plot_neutral_path(all_runs, jcss)
is_neutral = contains([all_runs.loading_condition], 'neutral', 'IgnoreCase', true);
neutral_path = all_runs(is_neutral);
specimen_names = unique([neutral_path.specimen]);
specimen_states = unique([neutral_path.state]);
colours = lines(numel(specimen_states));

for sn = 1:numel(specimen_names)
    % figure(sn)
    legend_text = "";

    current_specimen = neutral_path([neutral_path.specimen] == specimen_names(sn));

    for ss = 1:numel(specimen_states)
        current_state = current_specimen([current_specimen.state] == specimen_states(ss));
        opt_jcs = {current_state.(jcss)};

        colour = colours(ss, :);
        legend_text(end+1) = specimen_states(ss);
        for oj = 1:numel(opt_jcs)

            datum = opt_jcs{oj};
            flex = datum.Properties.VariableNames{strcmpi(datum.Properties.VariableNames, 'Flexion')};


            fieldnames = setdiff(datum.Properties.VariableNames, flex);
            for f = 1:numel(fieldnames)
                hold on;
                nexttile(f)
                fname = fieldnames{f};
                [x_arrowed, y_arrowed] = arrowed_line(datum.(flex), datum.(fname), 10, 100, 100);
                h = plot(x_arrowed, y_arrowed, 'Color', colour);
                xlabel("Flexion")
                ylabel(replace(fname, '_', ' '))
                grid on;

                if oj == 1 && f == 1
                    legend_handles(sn, ss) = h;
                end
            end
        end
    end
    sgtitle([specimen_names(sn) replace(jcss, '_', ' ')]);

    is_line = arrayfun(@(x) isa(x, 'matlab.graphics.chart.primitive.Line'), legend_handles(sn, :));
    current_legends = legend_handles(sn, :);
    legend_lines = current_legends(is_line);

    legend_text = legend_text(~(legend_text == ""));
    legend_text = state_regex_inv(legend_text);


    legend_text = legend_text(is_line);
    legend(legend_lines, legend_text)
end
end
