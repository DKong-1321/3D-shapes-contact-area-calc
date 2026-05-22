function visualise_tracker_transforms(gTf0, gTt0, is_right_knee)
f_x = gTf0(1:3,1); f_y = gTf0(1:3,2); f_z = gTf0(1:3,3);
    scatter3(f_x, f_y, f_z); hold on; grid
    t_x = gTt0(1:3,1); t_y = gTt0(1:3,2); t_z = gTt0(1:3,3);
    scatter3(t_x, t_y, t_z);
    xlabel("x"); ylabel("y"); zlabel("z");
    legend(["femur", "tibia"]);
    hold off;
% figure;
% hold on;
% grid;
% if is_right_knee
%     title("Right knee")
% else
%     title("Left knee")
% end
% %% Tibia
%     t = landmarks.tibia;
%     plot_t = scatter3(t.medial.Tx(1), t.medial.Ty(1), t.medial.Tz(1), 'b');
%     text(t.medial.Tx(1), t.medial.Ty(1), t.medial.Tz(1), "  Medial");
%     scatter3(t.lateral.Tx(1), t.lateral.Ty(1), t.lateral.Tz(1), 'b');
%     text(t.lateral.Tx(1), t.lateral.Ty(1), t.lateral.Tz(1), "  Lateral");
%     scatter3(t.distal.Tx(1), t.distal.Ty(1), t.distal.Tz(1), 'b');
%     text(t.distal.Tx(1), t.distal.Ty(1), t.distal.Tz(1), "  Distal");
% 
%     % Medial-lateral axis
%     o = (t.medial.translations_mean + t.lateral.translations_mean)/2;
%     scatter3(o(1), o(2), o(3), 8, 'b', 'filled')
%     if is_right_knee
%         med_lat = t.lateral.translations_mean - t.medial.translations_mean;
%         quiver3(t.medial.Tx(1), t.medial.Ty(1), t.medial.Tz(1), med_lat(1), med_lat(2), med_lat(3), 0, 'b')
%     else
%         med_lat = t.medial.translations_mean - t.lateral.translations_mean;
%         quiver3(t.lateral.Tx(1), t.lateral.Ty(1), t.lateral.Tz(1), med_lat(1), med_lat(2), med_lat(3), 0, 'b')
%     end
%     % Proximal-distal axis
%     prox_dist = t.distal.translations_mean - o;
%     quiver3(o(1), o(2), o(3), prox_dist(1), prox_dist(2), prox_dist(3), 0, 'b');
% 
% %% Femur
%     f = landmarks.femur;
% 
%     plot_f = scatter3(f.medial.Tx(1), f.medial.Ty(1), f.medial.Tz(1), 'r');
%     text(f.medial.Tx(1), f.medial.Ty(1), f.medial.Tz(1), '  Medial')
%     scatter3(f.lateral.Tx(1), f.lateral.Ty(1), f.lateral.Tz(1), 'r');
%     text(f.lateral.Tx(1), f.lateral.Ty(1), f.lateral.Tz(1), '  Lateral')
%     scatter3(f.proximal.Tx(1), f.proximal.Ty(1), f.proximal.Tz(1), 'r');
%     text(f.proximal.Tx(1), f.proximal.Ty(1), f.proximal.Tz(1), "  Proximal")
% 
%     % Medial-lateral axis
%     o = (f.medial.translations_mean + f.lateral.translations_mean)/2;
%     scatter3(o(1), o(2), o(3), 8, 'r', 'filled')
%     if is_right_knee
%         med_lat = f.lateral.translations_mean - f.medial.translations_mean;
%         quiver3(f.medial.Tx(1), f.medial.Ty(1), f.medial.Tz(1), med_lat(1), med_lat(2), med_lat(3), 0, 'r')
%     else
%         med_lat = f.medial.translations_mean - f.lateral.translations_mean;
%         quiver3(f.lateral.Tx(1), f.lateral.Ty(1), f.lateral.Tz(1), med_lat(1), med_lat(2), med_lat(3), 0, 'r')
%     end
%     % Proximal-distal axis
%     prox_dist = f.proximal.translations_mean - o;
%     quiver3(o(1), o(2), o(3), prox_dist(1), prox_dist(2), prox_dist(3), 0, 'r');
% 
% legend([plot_t, plot_f], {"Tibia", "Femur"})
% hold off;
end
