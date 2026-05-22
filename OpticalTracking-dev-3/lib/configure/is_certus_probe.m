function R = is_certus_probe(name)
% If there are multiple words, then it must be the probe. If not, then the word is one of the bone trackers.
% If you want to change how we determine that it's a probe, change this.
% Don't change the R = 'Probe' line.
is_probe = numel(strsplit(name, ' ')) > 1;

    if is_probe
        R = 'Probe';
    else
        R = name;
    end
end