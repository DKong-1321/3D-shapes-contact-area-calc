function trackers = polaris(data)
    probe_idx = find(contains(data.Properties.VariableNames, "Port"));

    n_probes = numel(probe_idx);
    n_tools = data.Tools(1);
    if n_tools == 1
        trackers = [];
        return
    end
    if n_probes == 1
        n_cols_between_probes = 0;
    else
        n_cols_between_probes = mode(diff(probe_idx));
    end
    names = data{1, probe_idx};
    % Exclude strays
    names = names(~contains(names, "Stray", "IgnoreCase", true));

    for k = 1:n_probes
        id = probe_idx(k);
        name = names{k};
        % Get the data and convert the errors into NaN
        datum = standardizeMissing(data(:, 1+id:id+n_cols_between_probes-1), -3.697314E+028);
        
        headers = datum.Properties.VariableNames;
        
        datum.Properties.VariableNames = clean_headers(headers);

        idx_tx_ty_tz = find(contains(headers, {'Tx', 'Ty', 'Tz'}), 3);

        q0 = datum.Q0;
        qx = datum.Qx;
        qy = datum.Qy;
        qz = datum.Qz;
        tx = datum{:, idx_tx_ty_tz(1)};
        ty = datum{:, idx_tx_ty_tz(2)};
        tz = datum{:, idx_tx_ty_tz(3)};
        error = datum.error;
        trackers(k) = Tracker(name, q0, qx, qy, qz, tx, ty, tz, error);
    end
end

function cleaned_headers = clean_headers(headers)
    headers = regexprep(headers, '.*Q0.*', 'Q0');
    headers = regexprep(headers, '.*Qx.*', 'Qx');
    headers = regexprep(headers, '.*Qy.*', 'Qy');
    headers = regexprep(headers, '.*Qz.*', 'Qz');
    headers = regexprep(headers, '.*Rx.*', 'Rx');
    headers = regexprep(headers, '.*Ry.*', 'Ry');
    headers = regexprep(headers, '.*Rz.*', 'Rz');
    headers = regexprep(headers, '.*Error.*', 'error');
    cleaned_headers = headers;
end
