function S = FEM_constant_line(local_coords, linear_elem, quadra_elem)
    % Linear
    S = zeros(3,2);
    
    %% Iterating at all nodes
    jacobian = linear_elem.jacobian(local_coords);
    p = 1;
    for gp1 = linear_elem.gauss
        w = gp1(2);

        % Linear Shape functions
        N1 = linear_elem.N(:,p)';

        % Quadratic Shape functions
        N2 = quadra_elem.N(:,p)';

        % Mass
        S = S + w * N2' * N1;
    end
    
    S  =  S * jacobian.det / 2;
end