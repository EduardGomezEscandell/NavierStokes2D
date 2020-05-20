function [M1, M12, M2, K, G1, G2] = FEM_constant(local_coords, linear_elem, quadra_elem)
    % Linear
    M1 = zeros(4,4);
    K = zeros(4,4);
    
    % Quadratic
    M2 = zeros(9,9);
    
    % Hybrid
    M12 = zeros(4,9);
    G1 = zeros(4,9);
    G2 = zeros(4,9);
    
    %% Iterating at all nodes
    
    jacobian1 = linear_elem.jacobian(local_coords);
    jacobian2 = quadra_elem.jacobian(local_coords);
    
    p = 1;
    for gp1 = linear_elem.gauss
        xi = gp1(1);
        for gp2 = linear_elem.gauss
            eta = gp2(1);
            w = gp1(2) * gp2(2);

            % Linear Shape functions
            J1 = jacobian1.calc(jacobian2,xi,eta);
            gradN1 = J1\linear_elem.gradN(:,:,p);
            N1 = linear_elem.N(:,p)';
            
            % Quadratic Shape functions
            J2 = jacobian2.calc(jacobian2,xi,eta);
            gradN2 = J2\quadra_elem.gradN(:,:,p);
            N2 = quadra_elem.N(:,p)';

            % Mass
            M1  =  M1  + w * (N1'*N1);
            M12 =  M12 + w * (N1'*N2);
            M2  =  M2  + w * (N2'*N2);
            
            % Diffusion
            K  =  K  + w * (gradN1' * gradN1);
            
            % G
            G1 = G1 - w * gradN1(1,:)' * N2;
            G2 = G2 - w * gradN1(2,:)' * N2;
        end
    end

    detJ = det(jacobian1.calc(jacobian1,0,0));
    M1  =  M1 * detJ;
    M12 = M12 * detJ;
    M2  =  M2 * detJ;
    K  =  K * detJ;
    G1 = G1 * detJ;
    G2 = G2 * detJ;
end