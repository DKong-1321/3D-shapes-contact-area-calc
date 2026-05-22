function self = add_transforms(self, label, transforms)
    self.Transform.(label) = transforms;
end
