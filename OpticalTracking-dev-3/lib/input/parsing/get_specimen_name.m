function specimen_name = get_specimen_name(root)
if contains(root, filesep)
    parent_folder_full_fp = split(root, filesep);
    test_name = parent_folder_full_fp{end};
else
    test_name = root;
end
    test_name_split= split(test_name, '_');
    specimen_name = test_name_split{1};
end