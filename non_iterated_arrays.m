function [M, K, K1, C1, G1, G2] = non_iterated_arrays(coords, X, visc, linear_elem, quadra_elem)  
    M = zeros(9,9);
    K = zeros(9,9);
    K1 = zeros(9,9);
    C1 = zeros(9,9);
    G1 = zeros(4,9);
    G2 = zeros(4,9);
    
    jacobian = quadra_elem.jacobian(coords);
            
    p = 1;
    for gp1 = quadra_elem.gauss
        xi = gp1(1);
        for gp2 = quadra_elem.gauss
            eta = gp2(1);
            w = gp1(2) * gp2(2);

            % Shape functions
            J = jacobian.calc(jacobian,xi,eta);
            N2 = quadra_elem.N(:,p)';
            gradN2 = J\quadra_elem.gradN(:,:,p);

            % Mass
            M  =  M  + w * (N2'*N2);

            % Diffusion
            K  =  K  + w * (gradN2' * gradN2);

            % Convection
            K1 = K1  + w * gradN2' * (N2*visc) * gradN2;
            C1 = C1  + w * N2' * ((N2*X(:,1:2) * gradN2));
        end
    end
    
    for gp1 = linear_elem.gauss
        xi = gp1(1);
        for gp2 = linear_elem.gauss
            eta = gp2(1);
            w = gp1(2) * gp2(2);

            % Shape functions
            N1 = linear_elem.N(:,p)';
            
            % G
            G1 = G1 + w * N1' * gradN2(1,:);
            G2 = G2 + w * N1' * gradN2(2,:);
        end
    end

    detJ = det(jacobian.calc(jacobian,0,0));
    M  =  M * detJ;
    K  =  K * detJ;
    K1 = K1 * detJ;
    C1 = C1 * detJ;
    G1 = G1 * detJ;
    G2 = G2 * detJ;
end