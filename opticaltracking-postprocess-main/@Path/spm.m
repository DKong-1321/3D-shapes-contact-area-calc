function spmi = spm(self)
    arguments
        self Path
    end
    flattened = self.flatten();
    spmi = SPM(flattened);
end
