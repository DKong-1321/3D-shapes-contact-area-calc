function root_ls = get_root_files(root, exclude_folders)
% `exclude_folders` is a cell with parts of names of folders that should be
% ignored. Typically "digit" (to exclude digitisation, digitise, digitize,
% etc) and "calibr" to exclude any variation of "calibration".

root_ls = dir(root);
root_ls = root_ls(~ismember({root_ls.name}, {'.', '..'}));
root_ls = root_ls(~contains({root_ls.name}, exclude_folders, 'IgnoreCase', true));
root_ls = Option(root_ls);
end