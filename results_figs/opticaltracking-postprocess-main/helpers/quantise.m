function [out, headers] = quantise(data)
    out = cell(numel(data), 1);
    headers = data{1}.Properties.VariableNames;
    idx_flexion = find(strcmpi(headers, 'flexion'));
    X = cellfun(@table2array, data, 'UniformOutput', false);
    X = cellfun(@(c) c(all(~isnan(c), 2), :), X, 'UniformOutput', false);

    % first = round(cellfun(@(x) x(1, idx_flexion), X));
    first = round(cellfun(@(x) min(x(:, idx_flexion)), X));
    last = round(cellfun(@(x) x(end, idx_flexion), X));
    [peak, peak_idx] = cellfun(@(x) max(x(:, idx_flexion)), X);
    peak = round(peak);

    flexion_arc = min(first):1:max(peak);
    if peak == last
        extension_arc = [];
    else
        extension_arc = max(peak):-1:min(last);
    end

    arc = [flexion_arc, extension_arc];

    all_groups = cellfun(@(x, pk_idx) compute_groups(x, pk_idx, idx_flexion, flexion_arc, peak), ...
                     X, num2cell(peak_idx), 'UniformOutput', false);

    n_rows = max([numel(arc) cellfun(@max, all_groups)]);
    for r = 1:numel(X)
        run = X{r};
        groups = all_groups{r};
        % quantised = accumarray(groups, (1:size(run,1))', [n_rows 1], @(x) {mean(run(x, :), 1)}, {nan(size(headers))});
        

        n_samples = size(run, 1);

        S = sparse(groups, 1:n_samples, 1, n_rows, n_samples);
        counts = full(sum(S, 2));

        quantised = (S * run) ./ counts;
        quantised(counts == 0, :) = NaN;
        
        % out{r} = array2table(vertcat(quantised{:}), "VariableNames", headers);
        out{r} = array2table(quantised, "VariableNames", headers);
    end
end
% function out = quantise(data)
%     headers = data{1}.Properties.VariableNames;
%     is_flexion = strcmpi(headers, 'flexion');
%     flexion = headers{is_flexion};
%
%     first = round(cellfun(@(x) x.(flexion)(1), data));
%     last = round(cellfun(@(x) x.(flexion)(end), data));
%     [peak, peak_idx] = cellfun(@(x) max(x.(flexion)), data);
%
%     peak = round(peak);
%     flexion_arc = min(first):1:max(peak);
%     extension_arc = max(peak):-1:min(last);
%
%
%     % idx_flex = find(strcmpi(headers, 'flexion'));
%     % X = cellfun(@table2array, data, 'UniformOutput', false);
%     for r = 1:numel(data)
%         run = data{r};
%         run_flex = run(1:peak_idx(r), :);
%         run_ext = run(peak_idx(r)+1:end, :);
%
%         out_flex = array2table(nan(numel(flexion_arc), size(run_flex, 2)), "VariableNames", headers);
%
%         out_ext = array2table(nan(numel(extension_arc), size(run_ext, 2)), "VariableNames", headers);
%         for n = 1:numel(flexion_arc)
%             ang = flexion_arc(n);
%             out_flex(n, :) = mean(run_flex(round(run_flex.(flexion)) == ang, :));
%         end
%         for n = 1:numel(extension_arc)
%             ang = extension_arc(n);
%             out_ext(n, :) = mean(run_ext(round(run_ext.(flexion)) == ang, :));
%         end
%
%         out{r} = [out_flex; out_ext];
%     end
% end

function groups = compute_groups(run, pk_idx, idx_flexion, flexion_arc, peak)
    flexion = round(run(:, idx_flexion));
    is_flexion = (1:size(run, 1))' <= pk_idx;
    groups = max(peak) - flexion + numel(flexion_arc);
    groups(1:pk_idx) = flexion(is_flexion) - min(flexion_arc) + 1;
end