function trackers = load_data(folder_path, config)
    csv_files = dir(fullfile(folder_path, "*.csv"));
    csv_files = csv_files(~contains({csv_files.name}, {'3d.csv', 'tension', 'sensor'})); % Certus data files that shouldn't be used
    tsv_files = dir(fullfile(folder_path, "*.tsv"));
    files = [csv_files; tsv_files];

    if isempty(files)
        trackers = Option.None;
        return
    end

    % Open the files as tables, use the headers to determine which optical
    % trackin system is being used, then apply any parsing necessary to
    % make data usable.
    % e.g. convert -3.697314E+028 => NaN on the polaris.
    for i = 1:numel(files)
        [~, parent_path, ~] = fileparts(files(i).name);
        fp = fullfile(files(i).folder, files(i).name);

        try
            data = readtable(fp, "VariableNamingRule","preserve");
        catch
            remove_whitespaces(fp);
            try
                data = readtable(fp, "VariableNamingRule","preserve");
            catch
                try
                    data = readtable(fp, "VariableNamingRule","preserve", "FileType","text", "Delimiter", '\t');
                catch
                    error("Malformed csv file: %s", fp)
                end
            end
        end


        tracker = Camera.load_data(data);
        if isempty(tracker)
            continue
        end
        trackers(:, i) = tracker;
        landmark = clean_specimen_condition(parent_path);
        trackers(:, i).add_landmark(landmark);
    end
    
    trackers.add_labels(config.camera_labels);

    trackers = Option.Some(trackers);
end

function remove_whitespaces(filepath)
%% Remove whitespaces, a common cause for matlab to not properly interpret the headers of a csv file.
% This is chatgpt stuff. Worth making sure it can't be done less painfully.
fid = fopen(filepath, 'r');
file_content = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);

file_content = file_content{1};

% Certus data always starts with Frame 1. We're looking for it to remove
% whitespaces. In the future, maybe just do it for every numeric row.
idx_start = 0;
for i = 1:10
    line = file_content{i};
    if startsWith(strtrim(line), '1')
        idx_start = i;
        break;
    end
end

% If it can't find a Frame 1 in the first 10 lines, then it's probably Polaris data and doesn't need to remove whitespaces
if idx_start == 0
    return;
end

file_content(i:end) = regexprep(file_content(i:end), '\s+', ''); % Remove all whitespace

% write the content back to a new file
fid = fopen(filepath, 'w');
if fid == -1
    error('Could not open the output file.');
end

for i = 1:length(file_content)
    fprintf(fid, '%s\n', file_content{i});
end

fclose(fid);
end
