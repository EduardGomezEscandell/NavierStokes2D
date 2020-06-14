function visc = get_viscosity(conc, visc0, derivative)
    switch derivative
        case 0
            %visc = visc0 * ones(size(conc));
            visc = visc0 + visc0 ./ (1 + exp(-10*(conc - 0.5)));
        case 1
            visc = 10 * visc0 * exp(-10*(conc - 0.5));
    end
end