function o = filter_state(self, state)
    error("Not yet implemented");
    mask = strcmpi(self.States, state);
    o = self(mask, :);
end
