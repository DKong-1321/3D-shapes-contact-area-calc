function probe = get_probe_labels(labels, calibration, root)
files = dir(fullfile(root, "*.tbr"));
if numel(files) == 0
    disp("x")
end
config = fullfile(files(1).folder, files(1).name);

probes = {calibration.probes};
names = {probes{1}.probe};
for i = 1:numel(names)
    probe(i).name = names{i};
    name = names{i}{1};
    label = read_config(config, name);
    % If it finds no label, attempt to determine it from probe name
    if ~isempty(label)
        probe(i).label = label;
        continue
    end

    labels_cell = struct2cell(labels);
    found_label_cells = cellfun(@(x) any(strcmp(x, strsplit(name, ' '))), labels_cell, 'UniformOutput',false);
    label_idx = cell2mat(found_label_cells);
    label = Option(labels_cell{label_idx});

    if label.is_none()
        non_empty = labels_cell(~cellfun(@isempty,labels_cell));
        idx = cellfun(@(x) any(contains(name, x)) , non_empty, 'UniformOutput',false);
        idx = [idx{:}];
        label = Option(non_empty(idx));
    end

    probe(i).label = label;
end
end

function label = read_config(file_path, target_string)
    % Open the file for reading
    label = [];
    fileID = fopen(file_path, 'r');

    if fileID == -1
        error('Failed to open the file.');
    end

    % Initialize variables
    found_line = '';
    string_found = false;

    % Read the file line by line and see if line contains target string 
    while ~feof(fileID)
        line = fgetl(fileID);

        % Check if the target string is at the start of the line
        if contains(line, target_string)
            probe = line;
            line = fgetl(fileID);

            found_line = strrep(line, target_string, '');  % Remove the target string
            string_found = true;
            break;  % Stop reading the file when the target string is found at the start of the line
        end
    end

    % Close the file
    fclose(fileID);

    % Check if the line was found or not
    if string_found
        % Split the found line into numeric values and convert into a 4x4 matrix
        values = regexp(found_line, '=', 'split');
        if contains(values, 'Label')
            fprintf('Found %s for probe %s\n', line, target_string);
            label = values(2:end);
        end
    end
end
