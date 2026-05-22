function full_runs = split_run(input, minima)
    headers = input.Properties.VariableNames;
    is_flexion = strcmpi(headers, 'flexion');
    flexion = headers{is_flexion};
    flex = input.(flexion);
    idx_minima = [1 minima length(flex)];

    n_possible_runs = numel(idx_minima) - 1;
    for n = 1:n_possible_runs
        start_idx = idx_minima(n);
        end_idx = idx_minima(n+1);
        runs{n} = input(start_idx:end_idx, :); 
    end

    includes_a_peak = cellfun(@(x) any(x.(flexion) > rms(flex, "omitnan")), runs);
    peak_flexion = cellfun(@(x) max(x.(flexion)), runs);
    is_complete_run = peak_flexion > 0.85*max(peak_flexion);
    full_runs = runs(is_complete_run & includes_a_peak);
    % 
    
    % 
    % run_with_flex_arc = runs(includes_a_peak);
end