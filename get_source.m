function source = get_source(velx, vely, derivative)
    switch derivative
        case 0
            %source = ones(size(vely));
            modU = sqrt(velx.^2 + vely.^2);
            source = 1./(1 + exp(-10*(modU - 0.5)));
        case 11 %ds/du
            modU = sqrt(velx.^2 + vely.^2);
            source = 1484.13 * velx .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
            source(isnan(source)) = 0;
        case 12 %ds/dv
            modU = sqrt(velx.^2 + vely.^2);
            source = 1484.13 * vely .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
            source(isnan(source)) = 0;
    end
end