classdef Dunnett
    properties
        signal
        dof
        state
        loading_condition
        p
        is_significant
    end
    methods
        function self = Dunnett(spm, control, p_critical, data)
            k = 1;
            if nargin > 0
            loading_conditions = spm.loading_conditions;
            states = spm.states;
            signal = string(spm.signal);
            dof = string(spm.dof);

            num_non_control_states = numel(states) - 1; % subtract control
            total_iter = numel(loading_conditions) * num_non_control_states;

            % Preallocate
            self(total_iter).loading_condition = [];

            
            for lc = 1:numel(loading_conditions)
                loading_condition = loading_conditions{lc};
                is_lc = loading_condition == [data.(signal).loading_condition];

                is_control_state = control == [data.(signal).state];
                is_control = is_control_state & is_lc;

                for st = 1:numel(states)
                    state = states{st};
                    if state == control
                        continue
                    end

                    is_state = state == [data.(signal).state];
                    mask = is_lc & is_state;
                    current = extract_trajectory_matrix(data.(signal)(mask), dof);
                    ctrl_spcm = extract_trajectory_matrix(data.(signal)(is_control), dof);

                    if all(current == 0 | isnan(current), "all")
                        continue
                    end

                    if size(ctrl_spcm, 2) ~= size(current, 2)
                        new_len = min(size(ctrl_spcm, 2), size(current, 2));
                        current = current(:, 1:new_len);
                        ctrl_spcm = ctrl_spcm(:, 1:new_len);
                    end

                    spm_t = spm1d.stats.ttest2(current, ctrl_spcm);
                    inference = spm_t.inference(p_critical);
                    if inference.h0reject
                        self(k).loading_condition = loading_condition;
                        self(k).state = state;
                        self(k).signal = signal;
                        self(k).p = inference.p;
                        self(k).is_significant = inference.z > inference.zstar;
                        self(k).dof = dof;
                        k = k + 1;
                    end
                end
            end
            end
            if k ~= 1
                self = self(1:k-1);
            end
        end

    end




end


function matrix = extract_trajectory_matrix(trajectories, dof)
    data = {trajectories.kinematics};
    data = cellfun(@(k) k.(dof)', data, 'UniformOutput', false);
    max_len = max(cellfun(@numel, data));
    matrix = cell2mat(cellfun(@(x) [x, nan(1, max_len - numel(x))], ...
                              data, 'UniformOutput', false)');
end