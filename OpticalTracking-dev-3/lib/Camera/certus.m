function trackers = certus(data)
    probe_idx = find(contains(data.Properties.VariableNames, "q0", "IgnoreCase",true));
    if isempty(probe_idx)
        error("Matlab was unable to read the headers from the csv. Try removing the spaces before the numbers.");
    end
    
    n = mode(diff(probe_idx)); % Columns 
    
    n_probes = length(probe_idx);
    labels = get_tracker_names(data.Properties.VariableNames(probe_idx));
    for k = 1:n_probes
        id = probe_idx(k);
    
        name = is_certus_probe(labels{k});
    
        datum = data(:, id:id+n-1);
        headers = datum.Properties.VariableNames;
        % Clean the headers, make them consistent across optical trackers
        datum.Properties.VariableNames = clean_headers(headers);

        q0 = datum.Q0;
        qx = datum.Qx;
        qy = datum.Qy;
        qz = datum.Qz;
        tx = datum.Tx;
        ty = datum.Ty;
        tz = datum.Tz;
        error = datum.Error;
        trackers(k) = Tracker(name, q0, qx, qy, qz, tx, ty, tz, error);
    end
end

function cleaned_headers = clean_headers(headers)
        headers = regexprep(headers, '.*q0.*', 'Q0');
        headers = regexprep(headers, '.*qx.*', 'Qx');
        headers = regexprep(headers, '.*qy.*', 'Qy');
        headers = regexprep(headers, '.*qz.*', 'Qz');
        headers = regexprep(headers, '.* x.*', 'Tx');
        headers = regexprep(headers, '.* y.*', 'Ty');
        headers = regexprep(headers, '.* z.*', 'Tz');
        cleaned_headers = regexprep(headers, '.*Error.*', 'Error');
end
    
function trackers = get_tracker_names(headers)
    split_names = cellfun(@(x) split(x, ' '), headers, 'UniformOutput', false);
    trackers = cellfun(@(x) strjoin(x(1:end-1), ' '), split_names, 'UniformOutput', false);
    trackers = lower(trackers);
end

function R = is_certus_probe(name)
% If there are multiple words, then it must be the probe. If not, then the word is one of the bone trackers.
% If you want to change how we determine that it's a probe, change this.
% Don't change the R = 'Probe' line.
is_probe = numel(strsplit(name, ' ')) > 1;

    if is_probe
        R = 'Probe';
    else
        R = name;
    end
end
