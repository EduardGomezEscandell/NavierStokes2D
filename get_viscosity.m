function visc = get_viscosity(conc, visc0)
    visc = ones(size(conc));
%     visc = visc0 + visc0 ./ (1 + exp(-10*(conc - 0.5)));
end