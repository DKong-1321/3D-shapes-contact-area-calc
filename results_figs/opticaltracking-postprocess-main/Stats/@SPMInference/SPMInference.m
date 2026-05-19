classdef SPMInference
    properties
        % rawdata
        signal
        dof
        within_subject
        between_subject
        states
        specimens
        loading_conditions
        alpha
    end
    methods
        function self = SPMInference(spm, alpha)
            if nargin > 0 
            n = numel(spm);
            self(n).signal = [];
            self(n).dof = [];
            self(n).within_subject = [];
            self(n).between_subject = [];
            self(n).states = [];
            self(n).specimens = [];
            self(n).loading_conditions = [];
            self(n).alpha = [];

            % self.rawdata = spm.rawdata;
            for i = 1:n
                self(i).signal = spm(i).signal;
                self(i).dof = spm(i).dof;
                self(i).states = spm(i).states;
                self(i).specimens = spm(i).specimens;
                self(i).alpha = alpha;
                self(i).loading_conditions = spm(i).loading_conditions;
                self(i).between_subject = spm(i).between_subject.inference(alpha);
                self(i).within_subject = spm(i).within_subject.inference(alpha);
            end
            end
        end
    end
end
