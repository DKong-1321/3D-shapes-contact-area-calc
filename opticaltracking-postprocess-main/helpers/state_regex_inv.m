function O = state_regex_inv(I)
    I = replace(I, '_w_', '+');
    I = replace(I, '_wo_', '-');
    O = replace(I, '_', ' ');
end
