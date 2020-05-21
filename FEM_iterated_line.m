function S = FEM_iterated_line(local_coords, visc, linear_elem, quadra_elem)
    % Linear
    S = zeros(3,3);
    
    %% Iterating at all nodes
    jacobian = linear_elem.jacobian(local_coords);
    p = 1;
    for gp1 = quadra_elem.gauss
        w = gp1(2);
        
        % Linear Shape functions
        N1 = linear_elem.N(:,p)';

        % Quadratic Shape functions
        N2 = quadra_elem.N(:,p)';
        invJ = jacobian.inv;
        gradN2 = invJ *quadra_elem.gradN(:,p)';

        % Mass
        S = S + w * N2' *  (N1*visc) * gradN2(1,:);
    end
    
    S  =  S * jacobian.det / 2;
end