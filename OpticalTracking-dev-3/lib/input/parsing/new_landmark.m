function landmark = new_landmark(trackers)
    landmark_names = {trackers.name};
    for i = 1:numel(trackers)
        landmark(i) = Landmark(landmark_names{i}, trackers);
    end
end
