function is_right_knee = get_knee_side(root, right_text, left_text)
    parent_folder_full_fp = split(root, filesep);
    test_name = parent_folder_full_fp{end};
    test_name_split = lower(string(split(test_name, {'_', ' ', '-'})));
    if any(test_name_split == lower(right_text), "all")
        is_right_knee = Option(true);
    elseif any(test_name_split == lower(left_text), "all")
        is_right_knee = Option(false);
    else
        files = string({dir(root).name});
        files = files(~ismember(lower(files), [".", "..", "results"]));
        % if isscalar(files)
        try
            is_right_knee = get_knee_side(files, right_text, left_text);
            is_right_knee = Option(is_right_knee);
        catch
            is_right_knee = Option.None;
        end
        % else
            % error("Could not determine knee side. Folder name should indicate the side.")
        % end
    end
end