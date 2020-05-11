function local_mat = FEM_matrices(coords, refelem, X, dX, visc)
    nodes_per_elem = length(coords);
    
    K = zeros(nodes_per_elem);
    K1 = zeros(nodes_per_elem);
    K21 = zeros(nodes_per_elem);
    K22 = zeros(nodes_per_elem);
    
    C1 = zeros(nodes_per_elem);
    C21 = zeros(nodes_per_elem);
    C22 = zeros(nodes_per_elem);
        
    M = zeros(nodes_per_elem);
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
                
                % Gradients
                grad_dRho= invJ * gradN * X(:,3);
                grad_du = invJ * gradN * dX(:,1);
                grad_dv = invJ * gradN * dX(:,2); 
                
                % Stiffnes
                K = K + w * (invJ*gradN)' * invJ*gradN;
                
                K1  = K1  + w * (invJ*gradN)' * (N*visc) * invJ*gradN;  % K1(nu)
                K21 = K21 + w * (invJ*gradN)' * grad_du * N;             % K2(u)
                K22 = K22 + w * (invJ*gradN)' * grad_dv * N;             % K2(v)
                
                C1  = C1  + w * N' * ((N*X(:,1:2) * invJ*gradN));  % C1(u,v)
                C21 = C21 + w * (N' * N) * grad_dRho(1);             % C21(rho)
                C21 = C21 + w * (N' * N) * grad_dRho(2);             % C22(rho)
                
                M = M + w * (N'*N);
                p = p+1;
            end
        end
        
        detJ = det(jacobian.calc(jacobian,0,0));
        local_mat.K   = K   * detJ;
        local_mat.K1  = K1  * detJ;
        local_mat.K21 = K21 * detJ;
        local_mat.K22 = K22 * detJ;
        local_mat.C1  = C1  * detJ;
        local_mat.C21 = C21 * detJ;
        local_mat.C22 = C22 * detJ;
        local_mat.M   = M   * detJ;
        
    else
       % Triangles 
    end
end