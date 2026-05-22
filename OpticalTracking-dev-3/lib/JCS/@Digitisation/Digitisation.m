classdef Digitisation < handle
    properties
        config
        trackers Tracker
        filepath
        bone
        transforms
        module
        angle_offset
    end
    methods (Static)
        function digitisations = new(root, config, module, angle) % => Digitisation
            arguments
                root
                config
                module Module
                angle = 0  % Digitisation angle. Leave as default if varying angles were used. It is possible to set maximum/minimum flexion angle on data instead
            end
            specimen_list = get_root_files(root, {'result', 'problem'}).unwrap(); % Get all files in root and exclude any folders that include `result`
            specimen_folders = fullfile({specimen_list.folder}, {specimen_list.name});
            for i = 1:numel(specimen_folders)
                % digitisations(i) = Option.None;
                specimen_name = get_specimen_name(specimen_list(i).name);
                config.specimen.name = specimen_name;
                fprintf("%d. Specimen: %s\n", i, specimen_name);
                filepath = get_root_files(specimen_folders{i}, {});
                if filepath.is_none
                    continue
                end
                filepath = filepath.unwrap();
                folder = get_folder(filepath, config.digitisation);
                root = folder.map(@fileparts);

                if root.is_none()
                    continue
                end
                root = root.unwrap();

                config.is_right_knee = get_knee_side(root, config.right, config.left) ...
                    .expect("Could not determine knee side. Check defaults.m => config.right and config.left for instructions");

                trackers = load_data(folder.unwrap(), config);
                if trackers.is_none()
                    continue
                end
                digitisation = Digitisation(trackers.unwrap(), root, config, module);
                switch module
                    case Module.Knee
                        Knee.assign_bone(digitisation);

                        [t, b] = Knee.calculate_transforms(digitisation.trackers(:, 1), "None", digitisation.transforms, config);
                        [digitisation_position.tf, digitisation_position.pf] = Knee.grood_and_suntay(b.femur, b.tibia, b.patella, t.fTt, t.fTp, config.is_right_knee);
                        signals = fields(digitisation_position);
                        for sg = 1:numel(signals)
                            signal = signals{sg};
                            if ~isempty(digitisation_position.(signal))
                                digitisation.angle_offset.(signal) = angle - mean(digitisation_position.(signal).flexion);
                            end
                        end
                    case Module.Hip
                        error("Not implemented")
                end

                if isempty(digitisation)
                    warning("Failed to complete digitisation");
                    continue
                end
                digitisations(i) = digitisation;

            end
        end


    end
    methods
        function path = root(self) % => [&str]
            paths = {self.filepath};
            path = fileparts(paths{1});
        end
    end

    methods (Access = private)
        function self = Digitisation(trackers, filepath, config, module) % => Digitisation
            self.trackers = trackers;
            self.config = config;
            self.filepath = filepath;
            self.module = module;
        end
    end
end

function path = get_folder(filepath, folders)
mask = contains({filepath.name}, folders , "IgnoreCase",true);
if sum(mask) > 1 % There's more than 1 digitisation file
    mask(1:find(mask, 1, "last")-1) = false; % Pick the last one
    warning("Found more than one digitisation folder. Using '%s'", filepath(mask).name);
end
path = fullfile(filepath(mask).folder, filepath(mask).name);
path = Option(path);
end
