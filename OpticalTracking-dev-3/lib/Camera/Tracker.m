classdef Tracker < handle
    properties
        Landmark
        Name
        Label
        Q0, Qx, Qy, Qz
        Rx, Ry, Rz,
        Tx, Ty, Tz,
        Error
    end
    properties (Access = private)
        Camera
    end
    methods (Static)
        function tracker = default(name, label)
            tracker = Tracker(name, [], [], [], [], [], [], []);
            tracker.Label = label;
        end
    end
    methods
        function self = Tracker(name, q0, qx, qy, qz, tx, ty, tz, error)
            self.Name  = name;

            self.Tx    = tx;
            self.Ty    = ty;
            self.Tz    = tz;

            self.Q0  = q0;
            self.Qx  = qx;
            self.Qy  = qy;
            self.Qz  = qz;

            [rx, ry, rz] = quaternion2euler(q0, qx, qy, qz);

            self.Rx    = rx;
            self.Ry    = ry;
            self.Rz    = rz;

            self.Error = error;
        end
        function camera = camera(self)
            camera = self.Camera;
        end
        function self = add_camera(self, camera)
            arguments
                self
                camera Camera
            end
            [self.Camera] = deal(camera);
        end
        function self = add_labels(self, labels)
            names = {self.Name};
            for i = 1:numel(names)
                name = names{i};
                labels_cell = self(i).Camera.get_possible_labels(labels);
                labels_cell = struct2cell(labels_cell);

                found_label_cells = cellfun(@(x) any(strcmp(x, strsplit(name, ' '))), labels_cell, 'UniformOutput',false);
                label_idx = cell2mat(found_label_cells);
                label = Option(labels_cell{label_idx});

                if label.is_none()
                    non_empty = labels_cell(~cellfun(@isempty,labels_cell));
                    idx = cellfun(@(x) any(contains(name, x)) , non_empty, 'UniformOutput',false);
                    idx = [idx{:}];
                    label = Option(non_empty(idx));
                end
                self(i).Label = label;
            end
        end

        function self = add_landmark(self, landmark)
            [self.Landmark] = deal(landmark);
        end
        function r = rotations(self)
            r = [self.Rx self.Ry self.Rz];
        end
        function r = rotations_mean(self)
            if isempty(self.rotations)
                r = [];
                return
            end
            r = mean(self.rotations(), "omitmissing");
        end
        function r = translations(self)
            r = [vertcat(self.Tx) vertcat(self.Ty) vertcat(self.Tz)];
        end
        function r = translations_mean(self)
            if isempty(self.translations)
                r = [];
                return
            end
            r = mean(self.translations(), "omitmissing");
        end
        function r = with_label(self, label)
            labels = [self.Label];
            has_label = labels.map(@(x) strcmp(x, label)).map(@cell2mat);
            if has_label.is_none()
                r = Option.None;
                return
            end
            labeled = self(has_label.unwrap());
            r = Option(labeled);
        end

        function result = contains(self, bone_position)
            % Looks for file names that include, e.g., "tibia" and "medial".
            % If none found, looks for TM (first letters capitalised).
            landmarks = {self.Landmark};

            % Match filename to whole word.
            bone = contains(landmarks, bone_position{1}, "IgnoreCase",true);
            pos = contains(landmarks, bone_position{2}, "IgnoreCase",true);
            if any(bone & pos)
                result = Option(self(bone & pos));
                return
            end

            one_letter = cellfun(@(x) x(1), bone_position);
            match_one_letter = strcmpi(landmarks, one_letter);
            if any(match_one_letter)
                result = Option(self(match_one_letter));
                return
            end

            % Attempt increasing lengths of the words
            substrings = @(word) arrayfun(@(n) word(1:n), length(word):-1:1, 'UniformOutput', false);
            all_substrings = cellfun(substrings, bone_position, 'UniformOutput', false);


            bone = contains(landmarks, all_substrings{1}, "IgnoreCase",true);
            pos = contains(landmarks, all_substrings{2}, "IgnoreCase",true);
            result = Option(self(bone & pos));
        end


    end
end



