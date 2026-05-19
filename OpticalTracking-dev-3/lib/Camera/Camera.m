classdef Camera
    enumeration
        Certus
        Polaris
        Unknown
    end
    methods (Static)
        function trackers = load_data(data)
            headers = data.Properties.VariableNames;
            camera = Camera.from_headers(headers);

            switch camera
                case Camera.Certus
                    trackers = certus(data);
                    if isempty(trackers)
                        trackers = [];
                        return
                    end
                    trackers.add_camera(Camera.Certus);
                case Camera.Polaris
                    trackers = polaris(data);
                    trackers.add_camera(Camera.Polaris);
                case Camera.Unknown
                    error("Unknown camera")
            end

            if numel(trackers) < 2
                warning("Fewer than two trackers found.");
            end
        end

        function camera = from_headers(headers)
            if strcmpi(headers{1}, 'frame')
                camera = Camera.Certus;
            elseif strcmpi(headers{1}, 'tools')
                camera = Camera.Polaris;
            else
                camera = Camera.Unknown;
            end
        end
    end

    methods
        function labels_cell = get_possible_labels(self, labels)
            switch self
                case Camera.Polaris
                    labels_cell = labels.polaris;
                case Camera.Certus
                    labels_cell = labels.certus;
            end
        end


    end
end
