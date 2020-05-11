function [M, K, K1, C1, Q1, Q2] = non_iterated_arrays(coords, X, visc, refelem)
    nodes_per_elem = length(coords);
    
    M = zeros(nodes_per_elem);
    K = zeros(nodes_per_elem);
    K1 = zeros(nodes_per_elem);
    C1 = zeros(nodes_per_elem);
    Q1 = zeros(nodes_per_elem);
    Q2 = zeros(nodes_per_elem);
    
    jacobian = refelem.jacobian(coords);
    
    if refelem.shape=='Q' % Quads
        
        % Integration
        p = 1;
        for gp1 = refelem.gauss
            xi = gp1(1);
            for gp2 = refelem.gauss
                eta = gp2(1);
                w = gp1(2) * gp2(2);
                
                % Shape functions
                J = jacobian.calc(jacobian,xi,eta);
                N = refelem.N(:,p)';
                gradN = J\refelem.gradN(:,:,p);
                
                % Mass
                M  =  M  + w * (N'*N);
                
                % Diffusion
                K  =  K  + w * (gradN' * gradN);
                
                % Convection
                K1 = K1  + w * gradN' * (N*visc) * gradN;
                C1 = C1  + w * N' * ((N*X(:,1:2) * gradN));
                Q1 = Q1 + w * N' * gradN(1,:);
                Q2 = Q2 + w * N' * gradN(2,:);
                
                p = p+1;
            end
        end
        
        detJ = det(jacobian.calc(jacobian,0,0));
        M  =  M * detJ;
        K  =  K * detJ;
        K1 = K1 * detJ;
        C1 = C1 * detJ;
        Q1 = Q1 * detJ;
        Q2 = Q2 * detJ;
    else
       % Triangles 
    end
end