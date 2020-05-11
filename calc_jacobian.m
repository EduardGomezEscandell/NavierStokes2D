function J = calc_jacobian(coords, visc0, theta, dt, X, K21, K22, C21, C22, Q1, Q2, T1, T2, M, A, B)
    n_nodes = size(coords,2);
    
    u_dof = 1:n_nodes;
    v_dof = n_nodes+1:2*n_nodes;
    p_dof = 2*n_nodes+1:3*n_nodes;
    d_dof = 3*n_nodes+1:4*n_nodes;

    
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    modU = modU + 1e-8; % Avoid errors in limit |U|->0
    
    % Derivtives of functions
    nu_prime = 10 * visc0 * exp(-10*(X(d_dof) - 0.5));
    s_of_u = 1484.13 * X(u_dof) .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
    s_of_v = 1484.13 * X(v_dof) .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
    
    % Turning into diagonal matrices
    Nr = spdiags(nu_prime,0,n_nodes,n_nodes);
    Su = spdiags(s_of_u,0,n_nodes,n_nodes);
    Sv = spdiags(s_of_v,0,n_nodes,n_nodes);
    
    % Jacobian
    Z = sparse(n_nodes, n_nodes);
    
    J = [               Q1                  Q2      Z                    Z
                         A                  Z      T1      dt*theta*K21*Nr
                         Z                  A      T2      dt*theta*K22*Nr
         dt*theta*C21-M*Su  dt*theta*C22-M*Sv       Z                   B];
     
end