function [K, C, M] = FEM_matrices(X, refelem, k,c,m)
    nodes_per_elem = length(X);
    
    K = zeros(nodes_per_elem, nodes_per_elem);
    C = zeros(nodes_per_elem, nodes_per_elem);
    M = zeros(nodes_per_elem, nodes_per_elem);
    jacobian = refelem.jacobian(X);
    
    if refelem.shape=='Q' % Quads
        
        % Integration
        p = 1;
        for xi = refelem.gauss_p
            for eta = refelem.gauss_p
                
                invJ = inv(jacobian.calc(jacobian,xi,eta));
                N = refelem.N(:,p)';
                gradN = refelem.gradN(:,:,p);
                
                K = K + (invJ*gradN)' * k * invJ*gradN;
                C = C + (invJ*gradN)' * c * N;
                M = M + m*(N'*N);
                
                p = p+1;
            end
        end
        
        area = det(jacobian.calc(jacobian,0,0));
        K = K * area/4;
        C = C * area/4;
        M = M * area/4;
        
    else
       % Triangles 
    end
end