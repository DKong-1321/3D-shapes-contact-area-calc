function path = path(self)
    names = unique([self.SpecimenName]);
    states = unique([self.SpecimenState]);
    loading_conditions = unique([self.LoadingCondition]);

    path = Path(self, names, states, loading_conditions);
    path.Root = self.Root;
end
