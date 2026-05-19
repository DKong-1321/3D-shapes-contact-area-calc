function path_normalised = normalise(self, neutral, native, fallback)
    arguments
        self Path
        neutral = "Neutral";
        native = "Native";
        fallback = "";
    end

    path_normalised = PathNormalised(self, neutral, native, fallback);
end
