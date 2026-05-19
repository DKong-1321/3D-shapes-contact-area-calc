function visualise_kinematics(transform)
N = length(transform.tibia);
fig = figure;
filename = 'tibia_rotation.gif'; % Output GIF file
for i = 1:N
    [tib] = unit_vectors(transform.tibia(:, :, i));
    [fem] = unit_vectors(transform.femur(:, :, i));
    clf(fig);
    grid on;
    zlabel('Flexion');
    ylabel('Varus/Valgus');
    xlabel('Internal/External');
    axis equal;
    title(["Rotations only. No correction for translations. Frame: " i]);
    hold on;
    %% Plot
    t = plotter(tib.x, tib.y, tib.z, tib.t, 'b');
    f = plotter(-fem.x, -fem.y, -fem.z, fem.t, 'r');
    legend([t f], {"Tibia", "Femur"})
    view(90+i/5,30);
    %%
    hold off
    % Capture the plot as a frame for the GIF
    frame = getframe(fig);
    % im = frame2im(frame);
    % [imind, cm] = rgb2ind(im, 256);
    % Write to GIF file
    % if i == 1
    %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    % else
    %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    % end
end
end

function rb = unit_vectors(transform)
if iscell(transform)
    T = transform{:};
else
    T = transform;
end
rb.x = T(1:3, 1);
rb.y = T(1:3, 2);
rb.z = T(1:3, 3);
rb.t = T(1:3, 4)/norm(T(1:3, 4));
end

function plot_x = plotter(x, y, z, t, colour)
origin = [0; 0; 0];
plot_x = quiver3(t(1,:), t(2,:), t(3,:), x(1,:), x(2,:), x(3,:), 0.2, colour);
quiver3(t(1,:), t(2,:), t(3,:), y(1,:), y(2,:), y(3,:), 0.2, colour);
quiver3(t(1,:), t(2,:), t(3,:), z(1,:), z(2,:), z(3,:), 0.2, colour); % Flexion
quiver3(t(1,:), t(2,:), t(3,:), x(1,:), y(2,:), z(3,:), 1, 'k');
axis([-1 1 -1 1 -2.5 0])
end
