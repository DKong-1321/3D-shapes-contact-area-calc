function digitisation = assign_bone(digitisation)
    % Uses file names and the config definitions (from defaults.m) to determine which file refers to which bone/location.
    % If you want to create this same setup for another set of bones, modify it as follows:
    % digitisation.bone.humerus.medial = trackers.contains({'humerus', 'medial'}).and_then(@(x) x.with_label(label.probe));
    %                   ^^^^^^^                             ^^^^^^^^^
    % Go into defaults.m and make sure there is a label for each of the new bones you want to introduce.
    arguments
        digitisation Digitisation
    end

    trackers = digitisation.trackers;
    config = digitisation.config;
    label = trackers.camera().get_possible_labels(config.camera_labels);
    %                                             ^^^^^^^^^^^^^^^^^^^^
    %                                        Looks for label definitions here

    digitisation.bone.tibia.medial = trackers.contains({'tibia', 'medial'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.tibia.lateral = trackers.contains({'tibia', 'lateral'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.tibia.distal = trackers.contains({'tibia', 'distal'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.tibia.tracker = trackers(:, 1).with_label(label.tibia);
    %% Femur
    digitisation.bone.femur.medial = trackers.contains({'femur', 'medial'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.femur.lateral = trackers.contains({'femur', 'lateral'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.femur.proximal = trackers.contains({'femur', 'proximal'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.femur.tracker = trackers(:, 1).with_label(label.femur);

    %% Patella
    digitisation.bone.patella.medial   = trackers.contains({'patella', 'medial'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.patella.lateral  = trackers.contains({'patella', 'lateral'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.patella.distal   = trackers.contains({'patella', 'distal'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.patella.inferior = trackers.contains({'patella', 'inferior'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.patella.superior = trackers.contains({'patella', 'superior'}).and_then(@(x) x.with_label(label.probe));
    digitisation.bone.patella.tracker  = trackers(:, 1).with_label(label.patella);

    % Calculate bone to tracker transforms in global coordinate system
    digitisation.transforms = Knee.bone_to_tracker_transform(digitisation.bone, config);
end
