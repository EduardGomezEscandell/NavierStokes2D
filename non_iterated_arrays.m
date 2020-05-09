function [K, K1, C1, M, b] = non_iterated_arrays(coords, X, dX, visc, refelem)
    nodes_per_elem = length(coords);
    
    K = zeros(nodes_per_elem);
    K1 = zeros(nodes_per_elem);
    C1 = zeros(nodes_per_elem);
    M = zeros(nodes_per_elem);
    b = zeros(nodes_per_elem);
    jacobian = refelem.jacobian(coords);
    
    if refelem.shape=='Q' % Quads
        
        % Integration
        p = 1;
        for xi = refelem.gauss_p
            for eta = refelem.gauss_p
                
                invJ = inv(jacobian.calc(jacobian,xi,eta));
                N = refelem.N(:,p)';
                gradN = refelem.gradN(:,:,p);
                w = 1; % placeholder
                
                K = K + w * (invJ*gradN)' * invJ*gradN;
                K1  = K1  + w * (invJ*gradN)' * (N*visc) * invJ*gradN;  % K1(nu)
                C1  = C1  + w * N' * ((N*X(:,1:2) * invJ*gradN));  % C1(u,v)
                M = M + w * (N'*N);
                p = p+1;
            end
        end
        
        detJ = det(jacobian.calc(jacobian,0,0));
        K   = K   * detJ;
        M   = M   * detJ;
        
    else
       % Triangles 
    end
end