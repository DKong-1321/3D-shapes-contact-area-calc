function p = plot(self)

    if isempty(self.Kinematics)
        p = plot(0);
        return
    end

    states = string(self.states);
    specimens = string(self.specimens);
    n_colours = length(fieldnames(self.Kinematics.(specimens(1))));
    colours = lines(n_colours);
    directions = string(self.directions);

    for sp = 1:numel(specimens)
        specimen = specimens(sp);
        for s = 1:numel(states)
            state = states(s);
            colour = colours(s, :);

            if any(cellfun(@(x) isempty(self.Kinematics.(specimen).(state).(x)), directions))
                continue
            end
            signals = fieldnames(self.Kinematics.(specimen).(state).(directions(1)));
            for sg = 1:numel(signals)
                figure(sg);
                signal = signals{sg};

                plots(s, sg) = gen_plots(self.Kinematics.(specimen).(state), directions, signal, colour);
            end
        end

        for sg = 1:numel(signals)
            signal = signals{sg};
            figure(sg);
            sgtitle(replace(signal, '_', ' '));
            legend(plots(:,sg), state_regex_inv(states));
        end
    end

end

function p = gen_plots(data, directions, jcs, colour)
    for d = 1:numel(directions)
        ap = directions(d);
        datum = data.(ap);
        if d > 1
            colour = 0.9 * colour;
        end

        if isempty(datum)
            p = plot(0);
            continue
        end



        orientations = datum.(jcs).Properties.VariableNames;
        orientations(contains(orientations, 'flexion')) = [];
        for o = 1:numel(orientations)
            nexttile(o); hold on;
            x = datum.(jcs).flexion;
            y = datum.(jcs).(orientations{o});
            [x_arrowed, y_arrowed] = arrowed_line(x, y, 10, 100, 100);
            % fill(x,y, colour, 'FaceAlpha', 0.1);
            p = plot(x_arrowed, y_arrowed, 'color', colour);
            % p = plot(x_arrowed, y_arrowed);
            grid on;
            axis square;
            xlabel("Flexion angle");
            ylabel(replace(orientations{o}, '_', ' '));
        end
    end
end
