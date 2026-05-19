function o = plot(self, orientations)
    arguments
        self
        orientations = [];
    end
    orientations = string(orientations);
    o = Plot(self, orientations);
end
% function o = plot(self, orientations)
%     if nargin > 1
%         orient = orientations;
%     else
%         orient = [];
%     end
%
%     if isempty(self.Kinematics)
%         o = plot(0);
%         return
%     end
%
%     states = self.States;
%     directions = self.Directions;
%     colours = lines(numel(states));
%     signals = self.Signals;
%
%     for sg = 1:numel(signals)
%         f(sg) = figure;
%         signal = signals{sg};
%         for s = 1:numel(states)
%             state = states(s);
%             colour = colours(s, :);
%
%             plots(s, sg) = gen_plots(self.Kinematics.(state).(signal), directions, colour, s, orient);
%
%         end
%         sgtitle(replace(signals(sg), '_', ' '));
%         legend(plots(:,sg), state_regex_inv(states));
%     end
% end
%
% function p = gen_plots(data, directions, colour, s, orientations)
%     arguments
%         data
%         directions
%         colour
%         s
%         orientations = [];
%     end
%     means = [];
%     for d = 1:numel(directions)
%         ap = directions(d);
%         datum = data.(ap);
%         if isempty(orientations)
%             orientations = datum.mean.Properties.VariableNames;
%         end
%         is_flexion = contains(orientations, 'flexion');
%         orientations(is_flexion) = [];
%         for o = 1:numel(orientations)
%             orientation = orientations{o};
%             means(o, d) = mean(datum.mean.(orientation));
%         end
%     end
%
%     is_first_higher = means(:, 1) > means(:, 2);
%
%     % linestyles = ["--", ":"];
%     for d = 1:numel(directions)
%         ap = directions(d);
%         datum = data.(ap);
%         if d > 1 %Differentiate anterior from posterior
%             colour = 0.9 * colour;
%         end
%
%         if isempty(datum.mean)
%             p = plot(0);
%             return
%         end
%
%
%
%         if isempty(orientations)
%             orientations = datum.mean.Properties.VariableNames;
%         end
%         is_flexion = contains(orientations, 'flexion');
%         orientations(is_flexion) = [];
%
%         step = 10;
%
%
%
%         for o = 1:numel(orientations)
%
%             nexttile(o); hold on;
%             x = datum.mean.flexion;
%             y = datum.mean.(orientations{o});
%             % p = plot(x, y, linestyles(d), 'color', colour);
%             p = plot(x, y, 'color', colour);
%
%             [orientation, multi] = invert_orientation(orientations{o});
%
%             y_std = multi * datum.std.(orientations{o});
%             idx = 1:step+1*s:numel(x);
%             is_bar_up = xor(is_first_higher(o), d > 1);
%             if is_bar_up
%                 errorbar(x(idx), y(idx), 0, y_std(idx), 'LineStyle', 'none', 'Color', colour*0.7);
%             else
%                 errorbar(x(idx), y(idx), y_std(idx), 0, 'LineStyle', 'none', 'Color', colour*0.7);
%             end
%             % fill(x,y, colour, 'FaceAlpha', 0.1);
%
%             grid on;
%             axis square;
%             xlabel("Flexion angle");
%             ylabel(replace(orientation, '_', ' '));
%         end
%
%
%     end
% end
%
% function [o, multi] = invert_orientation(orientation)
%     o = orientation;
%     multi = 1;
%     if strcmp(orientation, 'posterior')
%         o = 'anterior';
%         multi = -1;
%     end
% end
