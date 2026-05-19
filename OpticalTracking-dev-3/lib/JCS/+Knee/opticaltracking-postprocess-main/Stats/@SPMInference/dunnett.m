function out = dunnett(self, data, control)
    arguments
        self SPMInference
        data
        control = "Native"
    end
    if ~ismember(control, unique(vertcat(self.states)))
        error("Not a valid control. Options are: %s", strjoin(unique(vertcat(self.states)), ', '))
    end
    is_significant = @(c) cellfun(@(c) c.SPMs{3}.h0reject, c);
    between_subject = self(is_significant({self.between_subject}));
    within_subject = self(is_significant({self.within_subject}));


    for bs = 1:numel(between_subject)
        n_tests = numel(self(bs).states) - 1;

        p_critical = spm1d.util.p_corrected_bonf(self(bs).alpha, n_tests);


        dunn{bs} = Dunnett(between_subject(bs), control, p_critical, data.flatten());
        
    end
    out = [dunn{:}];
end
