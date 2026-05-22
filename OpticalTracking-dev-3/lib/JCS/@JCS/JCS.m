classdef JCS
    properties
        specimen
        loading_condition
        state
        transforms
        bones
        config
        module Module
        digitisation
    end
    properties(Access = private)
    end
    methods
        function jcs_solved = solve(self)
            % Alias for grood and suntay
            jcs_solved = [self.grood_and_suntay()];
        end
        function kinematics = grood_and_suntay(self)
            trajectories(numel(self)) = Trajectory();
            
            for i = 1:numel(self)
                switch self(i).module
                    case Module.Knee
                        right = self(i).config.is_right_knee;
                        fTt  = self(i).transforms.fTt;
                        fTp  = self(i).transforms.fTp;
                        femur = self(i).bones.femur;
                        tibia = self(i).bones.tibia;
                        patella = self(i).bones.patella;

                        trajectory = Trajectory(self(i).config.specimen.name, self(i).config.specimen.state, self(i).config.specimen.loading_condition, false, right);

                        [tf, pf] = Knee.grood_and_suntay(femur, tibia, patella, fTt, fTp, right);
                        % if ~isempty(tf)
                        %     tf.flexion = tf.flexion - self.digitisation.angle_offset.tf;
                        % end
                        % if ~isempty(pf)
                        %     pf.flexion = pf.flexion - self.digitisation.angle_offset.pf;
                        % end
                        trajectory.add_data('tibiofemoral', tf);
                        trajectory.add_data('patellofemoral', pf);

                        trajectory.add_transforms('gTfi', self(i).transforms.gTfi);
                        trajectory.add_transforms('gTti', self(i).transforms.gTti);
                        trajectory.add_transforms('gTpi', self(i).transforms.gTpi);
                        trajectory.add_transforms('fTt', self(i).transforms.fTt);
                        trajectory.add_transforms('fTp', self(i).transforms.fTp);

                        trajectory.set_root(self(i).digitisation.root);
                        
                        trajectories(i) = trajectory;
                    case Module.Hip
                        error("Not yet implemented");
                end
            end
            kinematics = Kinematics(trajectories, [self.digitisation], [self.config]);
        end
        function jcs_solved = helical_axis(self)
            error("Not yet implemented")
        end
        function jcs_solved = sara(self)
            error("Not yet implemented")
        end
    end

    methods (Static)
        function jcs = new(digitisations)
            i = 1;
            for d = 1:numel(digitisations)
                digitisation = digitisations(d);
                config = digitisation.config;
                digitisation_transforms = digitisation.transforms;
                module = digitisation.module;
                states = get_root_files(digitisation.filepath, config.digitisation).unwrap();
                for st = 1:numel(states)
                    state = states(st);
                    fp_data = fullfile(state.folder, state.name);
                    trackers = load_data(fp_data, config);
                    if trackers.is_none()
                        warning('No csv files in %s', fp_data)
                        continue;
                    end
                    trackers = trackers.unwrap();

                    config.specimen.state = clean_specimen_condition(state.name);

                    loading_conditions = unique({trackers.Landmark});
                    for lc = 1:numel(loading_conditions)
                        loading_condition = loading_conditions{lc};
                        config.specimen.loading_condition = loading_condition;
                        switch module
                            case Module.Knee
                                [transforms, bones] = Knee.calculate_transforms(trackers(:, lc), loading_condition, digitisation_transforms, config);
                            case Module.Hip
                                error("Not yet implemented");
                        end

                        jcs(i) = JCS(transforms, bones, digitisation, config);
                        jcs(i).module = module;
                        jcs(i).state = string(clean_specimen_condition(state.name));
                        jcs(i).specimen = string(config.specimen.name);
                        jcs(i).loading_condition = string(loading_condition);
                        i = i + 1;
                    end
                end
            end

            

        end
    end
    methods
        function self = JCS(transforms, bones, digitisation, config)
            self.bones = bones;
            self.transforms = transforms;
            self.config = config;
            self.digitisation = digitisation;
        end
    end
end
