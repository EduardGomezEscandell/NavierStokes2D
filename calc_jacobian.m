function J = calc_jacobian(coords, visc0, X, K21, K22, C21, C22, M, A, B)
    n_nodes = size(coords,2);
    
    u_dof = 1:n_nodes;
    v_dof = n_nodes+1:2*n_nodes;
    p_dof = 2*n_nodes+1:3*n_nodes;
    
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    modU = modU + 1e-8; % Avoid errors in limit |U|->0
    
    % Derivtives of functions
    nu_prime = 10 * visc0 * exp(-10*(X(p_dof) - 0.5));
    s_of_u = 1484.13 * X(u_dof) .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
    s_of_v = 1484.13 * X(v_dof) .* exp(10*modU) ./ (modU .* (148.413 + exp(10*modU)));
    
    % Turning into diagonal matrices
    Nr = spdiags(nu_prime,0,n_nodes,n_nodes);
    Su = spdiags(s_of_u,0,n_nodes,n_nodes);
    Sv = spdiags(s_of_v,0,n_nodes,n_nodes);
    
    % Jacobian
    empty_block = sparse(n_nodes, n_nodes);
    
    J = [               A        empty_block    0.5*K21*Nr;
              empty_block                  A    0.5*K22*Nr;
         0.5*C21-0.5*M*Su   0.5*C22-0.5*M*Sv             B];
     
end