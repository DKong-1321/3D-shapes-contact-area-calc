function o = filter_envelope(self, envelope)
    error("Not yet implemented");
    mask = contains(self.Directions, envelope, "IgnoreCase", true);
    if ~any(mask)
        o = [];
        return
    end

    o = self;
    o.Envelopes = self.Directions(mask);

    for s = 1:numel(o.States)
        state = o.States(s);
        env = fieldnames(o.Kinematics.(state));
        to_remove = setdiff(env, self.Directions(mask));
        o.Kinematics.(state) = rmfield(self.Kinematics.(state), to_remove);
    end
end
