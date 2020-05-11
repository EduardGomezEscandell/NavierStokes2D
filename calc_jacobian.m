function J = calc_jacobian(coords, connect, X, dX, visc0, mu, dt, removed_dof, refelem)
    n_nodes = size(coords,2);
    n_elems = size(connect,2);
    
    K_global = sparse(n_nodes,n_nodes);
    M_global = sparse(n_nodes,n_nodes);
    K1_global = sparse(n_nodes,n_nodes);
    K21_global = sparse(n_nodes,n_nodes);
    K22_global = sparse(n_nodes,n_nodes);
    C1_global = sparse(n_nodes,n_nodes);
    C21_global = sparse(n_nodes,n_nodes);
    C22_global = sparse(n_nodes,n_nodes);
    
    U = X(1:n_nodes);
    V = X(n_nodes+1:2*n_nodes);
    modU = sqrt(U.^2 + V.^2);
    
    rho = X(2*n_nodes+1:3*n_nodes);
    
    visc = visc0 + visc0 ./ (1 + exp(-10*(rho - 0.5)));
    
    for el=1:n_elems
        nodes = connect(:,el);
        local_coordinates = coords(:, nodes);

        x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
        dx = [dX(nodes),dX(nodes+n_nodes),dX(nodes+2*nodes)];
        
        v = visc(nodes,:);

        [K, K1, K21, K22, C1, C21, C22, M] = FEM_matrices(local_coordinates, refelem, x, dx, v);
        
        K_global(nodes,nodes) 	= K_global(nodes,nodes)   + K;
        M_global(nodes,nodes) 	= M_global(nodes,nodes)   + M;
        K1_global(nodes,nodes) 	= K1_global(nodes,nodes)  + K1;
        K21_global(nodes,nodes) = K21_global(nodes,nodes) + K21;
        K22_global(nodes,nodes) = K22_global(nodes,nodes) + K22;
        C1_global(nodes,nodes) 	= C1_global(nodes,nodes)  + C1;
        C21_global(nodes,nodes) = C21_global(nodes,nodes) + C21;
        C22_global(nodes,nodes) = C22_global(nodes,nodes) + C22;

    end
       
    % Shorter naming
    K = K_global;
    M = M_global;
    K1 = K1_global;
    K21 = K21_global;
    K22 = K22_global;
    C1 = C1_global;
    C21 = C21_global;
    C22 = C22_global;
    
    % Derivtives of functions
    nu_prime = 10 * visc0 * exp(-10*(rho - 0.5));
    s_of_u = 1484.13*U.*exp(10*modU) ./ (1e-8 + modU .* (148.413 + exp(10*modU)));
    s_of_v = 1484.13*V.*exp(10*modU) ./ (1e-8 + modU .* (148.413 + exp(10*modU)));
    
    Nr = diag(nu_prime);
    Su = diag(s_of_u);
    Sv = diag(s_of_v);
    
    empty_block = sparse(n_nodes, n_nodes);
    
    J = [ M / dt + 0.5*K1        empty_block             0.5*K21*Nr;
              empty_block    M / dt + 0.5*K1             0.5*K22*Nr;
         0.5*C21-0.5*M*Su   0.5*C22-0.5*M*Sv dt*M*0.5*mu*K + 0.5*C1];
     
end