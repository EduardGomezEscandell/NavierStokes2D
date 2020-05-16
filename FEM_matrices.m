function local_mat = FEM_matrices(coords, refelem, X, visc, mu, theta, dt)
    nodes_per_elem = 9;
    
    L = zeros(nodes_per_elem);
    supg = zeros(nodes_per_elem);
    
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
    
    detJ = det(jacobian.calc(jacobian,0,0));
    area = 4*detJ;
    h = sqrt(area);
    
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
                
                % Mass
                M = M + w * (N'*N);
                
                % Stabilization in concentration
                a = N*X(:,1:2);
                tau = (1/(theta*dt) + 2*norm(a)/h + 4*mu/h^2)^-1;
                supg = supg + tau * (a * gradN)' * N;
                
                % Stiffness
                K   =  K  + w * (gradN' * gradN);
                
                % Diffusion
                K1  =  K1 + w * gradN' * (N*visc) *gradN; % K1(nu)
                K21 = K21 + w * gradN' * grad_u * N;           % K2(u)
                K22 = K22 + w * gradN' * grad_v * N;           % K2(v)
                
                % Convection
                C1  =  C1 + w * N' * ((N*X(:,1:2) * gradN));  % C1(u,v)
                C21 = C21 + w * (N' * N) * grad_dRho(1);      % C21(rho)
                C22 = C22 + w * (N' * N) * grad_dRho(2);      % C22(rho)
                
                p = p+1;
            end
        end
        
        local_mat.supg= supg* detJ;
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