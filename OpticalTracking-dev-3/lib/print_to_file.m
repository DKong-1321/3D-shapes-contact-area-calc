function print_to_file(datum, path)
    tibiofemoral = setdiff(fieldnames(datum), "name");

    for i = 1:numel(tibiofemoral)
        fp = fullfile(path, "Results");
        mkdir(fp, tibiofemoral{i});
        for j = 1:numel(datum)
            loading_condition = datum(j).name;
            data = datum(j).(tibiofemoral{i});
            if ~isempty(data)
                filename = fullfile(fp, tibiofemoral{i}, [loading_condition '.csv']);
                writetable(data, filename);
            end
        end
    end
end
