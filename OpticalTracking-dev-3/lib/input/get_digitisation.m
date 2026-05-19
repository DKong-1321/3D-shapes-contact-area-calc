function fp_digitisation = get_digitisation(fp_conditions, digitisation)
    digitisation_file_mask = contains({fp_conditions.name}, digitisation , "IgnoreCase",true);
    if sum(digitisation_file_mask) > 1 % There's more than 1 digitisation file
        digitisation_file_mask(1:find(digitisation_file_mask, 1, "last")-1) = false; % Pick the last one
        warning("Found more than one digitisation folder. Using '%s'", fp_conditions(digitisation_file_mask).name);
    end
    fp_digitisation = fullfile(fp_conditions(digitisation_file_mask).folder, fp_conditions(digitisation_file_mask).name);
    fp_digitisation = Option(fp_digitisation);
end