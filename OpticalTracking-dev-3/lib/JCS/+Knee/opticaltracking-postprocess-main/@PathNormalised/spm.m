function spmlist = spm(self)
    arguments
        self PathNormalised
    end
    flattened = self.flatten();
    spmlist = [SPM(flattened)];
end
