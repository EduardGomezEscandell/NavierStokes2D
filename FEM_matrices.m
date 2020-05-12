function local_mat = FEM_matrices(coords, refelem, X, visc)
    nodes_per_elem = length(coords);    
    
    K = zeros(nodes_per_elem);
    K1 = zeros(nodes_per_elem);
    K21 = zeros(nodes_per_elem);
    K22 = zeros(nodes_per_elem);
    
    C1 = zeros(nodes_per_elem);
    C21 = zeros(nodes_per_elem);
    C22 = zeros(nodes_per_elem);
    
    Q1 = zeros(nodes_per_elem);
    Q2 = zeros(nodes_per_elem);
    
    M = zeros(nodes_per_elem);
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
                
                % Gradients
                grad_u = gradN * X(:,1);
                grad_v = gradN * X(:,2); 
                grad_dRho= gradN * X(:,4);
                
                % Stiffnes
                K   =  K  + w * (gradN' * gradN);
                
                K1  =  K1 + w * gradN' * (N*visc) *gradN; % K1(nu)
                K21 = K21 + w * gradN' * grad_u * N;           % K2(u)
                K22 = K22 + w * gradN' * grad_v * N;           % K2(v)
                
                C1  =  C1 + w * N' * ((N*X(:,1:2) * gradN));  % C1(u,v)
                C21 = C21 + w * (N' * N) * grad_dRho(1);           % C21(rho)
                C22 = C22 + w * (N' * N) * grad_dRho(2);           % C22(rho)
                
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