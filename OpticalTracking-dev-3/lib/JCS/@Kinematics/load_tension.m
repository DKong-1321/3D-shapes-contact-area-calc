function self = load_tension(self, identifier)
    arguments
        self
        identifier = "**/*tension.csv"
    end
    digitisation = [self.digitisation];
    path = digitisation.root;
    specimens = [self.trajectories.specimen];
    loading_conditions = [self.trajectories.LoadingCondition];
    states = [self.trajectories.states];


    files = Option(dir(fullfile(path, identifier)));

    if files.is_none
        output = Option.None();
        return
    end
    files = files.unwrap();

    for f = 1:numel(files)
        filename = fullfile(files(f).folder, files(f).name);

        [~, this_loading_condition, ~] = fileparts(files(f).name);
        this_loading_condition = clean_specimen_condition(this_loading_condition);
        this_loading_condition = replace(this_loading_condition, '_tension', '');
        words = split(filename, filesep);
        this_state = words{end-1};
        this_specimen_full = words{end-2};
        this_specimen = get_specimen_name(this_specimen_full);


        is_specimen = this_specimen == specimens;
        is_lc = this_loading_condition == loading_conditions;
        is_state = this_state == states;



        tension_raw = readtable(filename, "VariableNamingRule","preserve");

        % Get the force and z (flexion)
        tension = table();
        z = unwrap_windowed(tension_raw.Z);
        
        z(z > 400) = NaN; % Sometimes unwrap thinks it's 2 phases away, so this just removes them.
        if diff(z) > 100
            keyboard
        end
        tension.flexion = smoothdata(z, "gaussian", 8); % Fairly low number, but found it to copy the max/min best.
        forces = tension_raw.("Force (N)");
        if isnumeric(forces)
            tension.tension = forces;
        else
            tension.tension = str2double(erase(forces, ' N'));
        end
        
        % if any(diff(tension.force) > 30)
        %     keyboard
        % end

        mask = is_specimen & is_lc & is_state;
        if ~any(mask)
            continue
        end
        if sum(mask) > 1
            idx = find(mask);
            for i = 1:numel(idx)
                self.trajectories(idx(i)).add_sensor('tension', tension);
            end
        else
        self.trajectories(mask).add_sensor('tension', tension);
        end
   
    end
end

function y = unwrap_windowed(x, window_size, threshold)
    arguments
        x
        window_size = 5;
        threshold = 150;
    end

    x = x(:);

    y = zeros(size(x));
    y(1) = x(1);

     for k = 2:length(x)
        idx_start = max(1, k - window_size);
        ref = median(y(idx_start:k-1));

        delta = x(k) - ref;

        if delta > threshold
            delta = delta - 360;
        elseif delta < -threshold
            delta = delta + 360;
        end
        y(k) = ref + delta;
    end

end
