function probe = probe_labels_from_name(labels, calibration)
    names = unique({calibration.Name});
    for i = 1:numel(names)
        name = names{i};
        probe(i).name = name;
    
        labels_cell = struct2cell(labels);
        found_label_cells = cellfun(@(x) any(strcmp(x, strsplit(name, ' '))), labels_cell, 'UniformOutput',false);
        label_idx = cell2mat(found_label_cells);
        label = Option(labels_cell{label_idx});
    
        if label.is_none()
            non_empty = labels_cell(~cellfun(@isempty,labels_cell));
            idx = cellfun(@(x) any(contains(name, x)) , non_empty, 'UniformOutput',false);
            idx = [idx{:}];
            label = Option(non_empty(idx));
        end
    
        probe(i).label = label;
    end
end