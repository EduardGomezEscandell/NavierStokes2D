function source = get_source(velx, vely)
%       source = ones(size(vely));
    modU = sqrt(velx.^2 + vely.^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
end