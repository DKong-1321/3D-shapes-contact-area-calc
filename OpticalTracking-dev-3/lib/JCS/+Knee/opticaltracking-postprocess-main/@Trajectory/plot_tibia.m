function plot_tibia(self)
    range = 1:10:90;
    scale = 1;
    is_neutral = contains([self.LoadingCondition], "neutral", "IgnoreCase", true);
    is_native = contains([self.SpecimenState], ["native", "uka_w_acl"], "IgnoreCase", true);
    is_native_neutral = is_neutral & is_native;

    filtered = self(is_native_neutral);
    specimens = unique([filtered.SpecimenName]);
    states = unique([filtered.SpecimenState]);

    colours = lines(numel(states));
    for sp = 1:numel(specimens)
        figure; hold on;
        specimen = specimens(sp);
        for st = 1:numel(states)
            state = states(st);
            mask = [filtered.SpecimenName] == specimen & [filtered.SpecimenState] == state;
            transforms = filtered(mask).Transform;
            for t = 1:numel(transforms)
                headers = fields(transforms(t));
                for h = 1:numel(headers)
                    header = headers{h};

                    tTf = transforms.(header);
                    origins = squeeze(tTf(1:2, 4, :));
                    fTfx = [1; 0; 0; 1]; % Epicondylar axis in femur == x axis.
                    tTfx = squeeze(pagemtimes(tTf, fTfx)); %Epicondylar axis in tibial frame of reference;

                    x0 = origins(1, range);
                    y0 = origins(2, range);
                    u = tTfx(1, range);
                    v = tTfx(2, range);

                    x1 = x0 + scale * u;
                    y1 = y0 + scale * v;

                    X = [x0; x1];
                    Y = [y0; y1];

                    plot(X, Y, 'Color', colours(st, :), 'LineWidth',1.5);
                    text(x1, y1, string(range-1));
                    xlabel('x'); ylabel('y');
                    title([specimen; 'Epicondylar axis projected on tibial plateau'])
                end
            end
        end
    end
end
