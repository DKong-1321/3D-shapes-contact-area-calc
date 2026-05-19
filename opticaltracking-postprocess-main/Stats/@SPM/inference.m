function spmi =  inference(self, alpha)
    arguments
        self SPM
        alpha = 0.05
    end
    spmi = SPMInference(self, alpha);
end