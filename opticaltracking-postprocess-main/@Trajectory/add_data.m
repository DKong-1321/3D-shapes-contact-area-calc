function self = add_data(self, label, data)
    headers = data.Properties.VariableNames;
    if ~any(contains(headers, 'flexion'))
        error("Kinematics table must contain field 'flexion'.")
    end
    self.Kinematics.(label) = data;
end
