classdef Plot
    properties
        LineHandles
        States
        Signals
        Orientations
        FigureHandles
    end

    methods
        function self = Plot(envelope,orientations)
            arguments
                envelope EnvelopeAverage
                orientations = [];
            end
            self.States = envelope.States;
            self.Signals = envelope.Signals;
            self.Orientations = envelope.Directions;
            [self.LineHandles, self.FigureHandles] = create_plot(envelope, orientations);
        end
    end
end

function [line_handles, figure_handles] = create_plot(self, orientations)
    if isempty(self.Kinematics)
        line_handles = plot(0);
        return
    end

    states = self.States;
    directions = self.Directions; %anterior-posterior
    colours = lines(numel(states));
    signals = self.Signals;
    %orientations = degree of freedom

    for sg = 1:numel(signals)
        figure_handles(sg) = figure;
        signal = signals{sg};
        for s = 1:numel(states)
            state = states(s);
            colour = colours(s, :);

            plots(s).(signal) = gen_plots(self.Kinematics.(state).(signal), directions, colour, s, orientations);

        end
        sgtitle(replace(signals(sg), '_', ' '));
        plot_handles = [plots.(signal)];
        orient = fields(plot_handles);
        h = [plot_handles.(orient{1})];
        legend([h.(directions(1))], state_regex_inv(states));
    end
    line_handles = plots;
    
end

function p = gen_plots(data, directions, colour, s, DOFs)
    arguments
        data
        directions % directions of the envelope. "ant" and "pos", etc
        colour
        s
        DOFs = []; % posterior, medial, internal, etc
    end
    means = [];
    for d = 1:numel(directions)
        direction = directions(d);
        try
        datum = data.(direction);
        catch ME
            keyboard
        end
        if isempty(DOFs)
            DOFs = datum.mean.Properties.VariableNames;
        end
        is_flexion = contains(DOFs, 'flexion');
        DOFs(is_flexion) = [];
        for o = 1:numel(DOFs)
            dof = DOFs{o};
            means(o, d) = mean(datum.mean.(dof));
        end
    end

    is_first_higher = means(:, 1) > means(:, 2);

    % linestyles = ["--", ":"];
    for d = 1:numel(directions)
        direction = directions(d);
        datum = data.(direction);
        if d > 1 %Differentiate anterior from posterior
            colour = 0.9 * colour;
        end

        if isempty(datum.mean)
            for o = 1:numel(DOFs)
                p.(DOFs{o}) = plot(0);
            end
            return
        end



        if isempty(DOFs)
            DOFs = datum.mean.Properties.VariableNames;
        end
        is_flexion = contains(DOFs, 'flexion');
        DOFs(is_flexion) = [];

        step = 10;



        for o = 1:numel(DOFs)
            dof = DOFs{o};

            nexttile(o); hold on;
            x = datum.mean.flexion;
            y = datum.mean.(dof);
            % p = plot(x, y, linestyles(d), 'color', colour);
            p.(dof).(direction) = plot(x, y, 'color', colour);

            y_std = datum.std.(dof);
            idx = 1:step+1*s:numel(x);
            is_bar_up = xor(is_first_higher(o), d > 1);
            if is_bar_up
                errorbar(x(idx), y(idx), 0, y_std(idx), 'LineStyle', 'none', 'Color', colour*0.7);
            else
                errorbar(x(idx), y(idx), y_std(idx), 0, 'LineStyle', 'none', 'Color', colour*0.7);
            end
            % fill(x,y, colour, 'FaceAlpha', 0.1);

            grid on;
            axis square;
            xlabel("Flexion angle");
            ylabel(replace(dof, '_', ' '));
        end


    end
end
