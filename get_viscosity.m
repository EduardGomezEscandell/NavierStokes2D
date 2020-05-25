function visc = get_viscosity(conc, visc0)
%     visc = visc0 * ones(size(conc));
    visc = visc0 * (1 + conc);
    visc = visc0 + visc0 ./ (1 + exp(-10*(conc - 0.5)));
end