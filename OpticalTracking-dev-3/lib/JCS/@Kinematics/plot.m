function fig = plot(self, idx_interpol)
    arguments
    self JCS
    idx_interpol = [];
    end
    fig = self.trajectories.plot(idx_interpol);
end
