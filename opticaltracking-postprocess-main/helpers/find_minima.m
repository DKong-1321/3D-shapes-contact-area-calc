function minima = find_minima(input)
    if isempty(input)
        return
    end
    local_minima = islocalmin(input);
    if ~any(local_minima)
        minima = [];
        return
    end

    distance_threshold = round(10+7*log(height(input)));  % Minimum distance between minima


    % We get several minima on a long tail where the user forgot to turn off
    % the camera. We require a minimum slope to make sure it's actual data.
    idx_minima = find(local_minima);
    idx_sweep = idx_minima * [1 1 1] + [-round(.8*distance_threshold) 0 round(.8*distance_threshold)];
    idx_sweep(idx_sweep > numel(local_minima)) = idx_minima(end);
    idx_sweep(idx_sweep < 1) = 1;
    changes_over_time = diff(input(idx_sweep), 2, 2) > 1;

    large_change = false(size(local_minima));
    large_change(idx_minima(changes_over_time)) = true;

    close_to_extension = input < mean(input, "omitmissing");
    reasonable_minima = local_minima & large_change & close_to_extension;

    idx_reasonable_minima = find(reasonable_minima);
    % Find the minima that aren't too close to each other:

    distances = diff(idx_reasonable_minima);

    % With noisy data you end up with a few local minima. We do 2 passes to
    % find true minima: first one to find points that are far away from each
    % other, and second to pick a mid point between the leftover close ones.
    close_groups = logical([0; distances < distance_threshold]);
    if ~any(close_groups)
        minima = idx_reasonable_minima;
    else
        minima_clusters = {};          % Cell array to hold groups
        cluster = [];   % Temporary holder for the current group
        for j = 2:length(close_groups)

            % Create clusters of minima
            current = close_groups(j);
            previous = close_groups(j-1);
            if current == 1 && previous == 0
                cluster = union(idx_reasonable_minima(j-1), idx_reasonable_minima(j));
            elseif current == 1
                cluster = union(cluster, idx_reasonable_minima(j));
                % Flush cluster of minima when finding a new cluster.
            elseif current == 0
                if isempty(cluster)
                    cluster = idx_reasonable_minima(j-1);
                    minima_clusters{end+1} = cluster;
                end
                %start new cluster; flush the previous
                minima_clusters{end+1} = cluster;
                cluster = idx_reasonable_minima(j);
            end
        end

        % Missed the last cluster for some reason
        if current
            last_one = union(cluster, idx_reasonable_minima(j));
        else
            last_one = idx_reasonable_minima(j);
        end
        if ~isempty(minima_clusters)
            if isempty(minima_clusters) || isempty(intersect(minima_clusters{end}, last_one))
                minima_clusters{end+1} = last_one;
            end
        end

        minima = get_min_from_clusters(minima_clusters, input);
    end

    minima = nearby_if_nan(minima, input);
    minima = minima(:).';
end

function output = get_min_from_clusters(cluster, input)
output = [];
for i = 1:length(cluster)
    % Find the median value in the cluster of minima. Alternative is to find
    % the minimum on the cluster.
    med = median(cluster{i});
    output(i) = round(med);
    % [min_val, min_idx] = min(input(cluster{i}));
    % output(i) = cluster{i}(min_idx);
end
end

function minima = nearby_if_nan(minima, input)
% Moves left and right a few times until it finds something that isn't NaN.
if ~any(isnan(minima))
    return
end
for i = 1:length(minima)
    if isnan(input(minima(i)))
        left = minima(i) - 1;
        right = minima(i) + 1;

        while (left >= 1 && isnan(input(left))) && (right <= length(input) && isnan(input(right)))
            left = left - 1;
            right = right + 1;
        end

        if left >= 1 && ~isnan(input(left))
            minima(i) = left;
        elseif right <= length(input) && ~isnan(input(right))
            minima(i) = right;
        end
    end
end
end
