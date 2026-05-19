function visualise_landmarks_fixed(bone, config)
% Run directly: visualise_landmarks_fixed(digitisation(1).bone, digitisation(1).config)

figure('Name', 'Digitisation Landmarks', 'Color', 'w');
hold on; grid on; axis equal; view(3);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

%% Tibia (blue)
t     = bone.tibia;
t_med  = t.medial.unwrap().translations_mean(:);
t_lat  = t.lateral.unwrap().translations_mean(:);
t_dist = t.distal.unwrap().translations_mean(:);

scatter3(t_med(1),  t_med(2),  t_med(3),  80, 'b', 'filled');
text(t_med(1),  t_med(2),  t_med(3),  '  T-Medial',  'Color','b');
scatter3(t_lat(1),  t_lat(2),  t_lat(3),  80, 'b', 'filled');
text(t_lat(1),  t_lat(2),  t_lat(3),  '  T-Lateral', 'Color','b');
scatter3(t_dist(1), t_dist(2), t_dist(3), 80, 'b', 'filled');
text(t_dist(1), t_dist(2), t_dist(3), '  T-Distal',  'Color','b');

t_origin = (t_med + t_lat) / 2;
scatter3(t_origin(1), t_origin(2), t_origin(3), 60, 'b', 'filled');

if config.is_right_knee
    t_ml = t_lat - t_med;
    quiver3(t_med(1), t_med(2), t_med(3), t_ml(1), t_ml(2), t_ml(3), 0, 'b', 'LineWidth', 2);
else
    t_ml = t_med - t_lat;
    quiver3(t_lat(1), t_lat(2), t_lat(3), t_ml(1), t_ml(2), t_ml(3), 0, 'b', 'LineWidth', 2);
end
t_pd = t_origin - t_dist;
quiver3(t_dist(1), t_dist(2), t_dist(3), t_pd(1), t_pd(2), t_pd(3), 0, 'b', 'LineWidth', 2);

%% Femur (red)
f      = bone.femur;
f_med  = f.medial.unwrap().translations_mean(:);
f_lat  = f.lateral.unwrap().translations_mean(:);
f_prox = f.proximal.unwrap().translations_mean(:);

scatter3(f_med(1),  f_med(2),  f_med(3),  80, 'r', 'filled');
text(f_med(1),  f_med(2),  f_med(3),  '  F-Medial',   'Color','r');
scatter3(f_lat(1),  f_lat(2),  f_lat(3),  80, 'r', 'filled');
text(f_lat(1),  f_lat(2),  f_lat(3),  '  F-Lateral',  'Color','r');
scatter3(f_prox(1), f_prox(2), f_prox(3), 80, 'r', 'filled');
text(f_prox(1), f_prox(2), f_prox(3), '  F-Proximal', 'Color','r');

f_origin = (f_med + f_lat) / 2;
scatter3(f_origin(1), f_origin(2), f_origin(3), 60, 'r', 'filled');

if config.is_right_knee
    f_ml = f_lat - f_med;
    quiver3(f_med(1), f_med(2), f_med(3), f_ml(1), f_ml(2), f_ml(3), 0, 'r', 'LineWidth', 2);
else
    f_ml = f_med - f_lat;
    quiver3(f_lat(1), f_lat(2), f_lat(3), f_ml(1), f_ml(2), f_ml(3), 0, 'r', 'LineWidth', 2);
end
f_pd = f_prox - f_origin;
quiver3(f_origin(1), f_origin(2), f_origin(3), f_pd(1), f_pd(2), f_pd(3), 0, 'r', 'LineWidth', 2);

%% Title and legend
if config.is_right_knee
    title([config.specimen.name ' — Right knee — Digitisation landmarks']);
else
    title([config.specimen.name ' — Left knee — Digitisation landmarks']);
end
pt = scatter3(nan,nan,nan,60,'b','filled');
pf = scatter3(nan,nan,nan,60,'r','filled');
legend([pt pf], {'Tibia','Femur'});
hold off;
end