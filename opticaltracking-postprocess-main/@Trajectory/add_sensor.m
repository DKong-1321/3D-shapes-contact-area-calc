function self = add_sensor(self, label, data)
    headers = data.Properties.VariableNames;
    if ~any(contains(headers, 'flexion'))
        error("Sensor table must contain field 'flexion'.")
    end
    self.Sensors.(label) = data;
end
