function [K, K_v, C, M, F, Fprev] = FEM_matrices(X, refelem, a, s, s_prev, visc)
    nodes_per_elem = length(X);
    
    K = zeros(nodes_per_elem);
    K_v = zeros(nodes_per_elem);
    C = zeros(nodes_per_elem);
    M = zeros(nodes_per_elem);
    F = zeros(nodes_per_elem, 1 );
    Fprev = zeros(nodes_per_elem, 1 );
    jacobian = refelem.jacobian(X);
    
    if refelem.shape=='Q' % Quads
        
        % Integration
        p = 1;
        for xi = refelem.gauss_p
            for eta = refelem.gauss_p
                
                invJ = inv(jacobian.calc(jacobian,xi,eta));
                N = refelem.N(:,p)';
                gradN = refelem.gradN(:,:,p);
                w = 1; % placeholder
                
                K_v = K_v + w * (invJ*gradN)' * (N*visc) * invJ*gradN;
                K = K + w * (invJ*gradN)' * invJ*gradN;
                C = C + w * N' * ((N*a) * invJ*gradN);
                M = M + w * (N'*N);
                F = F + w * (N*s) * N';
                Fprev = Fprev + w * (N*s_prev) * N';
                p = p+1;
            end
        end
        
        detJ = det(jacobian.calc(jacobian,0,0));
        K = K * detJ;
        C = C * detJ;
        M = M * detJ;
        
    else
       % Triangles 
    end
end