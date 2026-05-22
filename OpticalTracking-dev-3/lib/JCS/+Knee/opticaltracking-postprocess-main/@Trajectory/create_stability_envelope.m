function envelope = create_stability_envelope(self, envelopes, native, neutral)
    arguments
        self
        envelopes
        native = "Native"
        neutral = "Neutral"
    end

    % names = unique([self.SpecimenName]);
    states = unique([self.SpecimenState]);

    is_native = contains([self.SpecimenState], native, "IgnoreCase", true);
    native_state = self(is_native);
    if isempty(native_state)
        error("Could not find native/intact state under the name %s. Try one of these: %s", native, strjoin(unique(self.states), ', '))
    end
    is_passive_flex = contains([native_state.LoadingCondition], neutral, "IgnoreCase", true);
    native_passive_flex = native_state(is_passive_flex);

    if isempty(native_passive_flex)
        error("Native Neutral flexion was not detected")
    end

    envelope = Envelope(self, envelopes, native_passive_flex, states);
end
