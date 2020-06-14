function J = build_jacobian(mesh, dof, bc_data, X, dt, theta, visc0, Au, Av, B, M12, K21, K22, C21, C22, G1, G2, T1, T2)   
    %% Derivatives of functions
    nu_prime = get_viscosity(X(dof.d), visc0, 1);
    s_of_u = get_source(X(dof.u), X(dof.v), 11);
    s_of_v = get_source(X(dof.u), X(dof.v), 12);
    
    % Turning into diagonal matrices
    Nr = spdiags(nu_prime,0,mesh.corners,mesh.corners);
    Su = spdiags(s_of_u,0,mesh.nodes,mesh.nodes);
    Sv = spdiags(s_of_v,0,mesh.nodes,mesh.nodes);    
    
    % Reducing
    Nr = Nr(bc_data.unkn_d, bc_data.unkn_d);
    
    Su  = Su(bc_data.unkn_u, bc_data.unkn_u);
    Sv  = Sv(bc_data.unkn_v, bc_data.unkn_v);
    
    %% Reducing matrices not reduced before
    K21 = K21(bc_data.unkn_u, bc_data.unkn_d);
    K22 = K22(bc_data.unkn_v, bc_data.unkn_d);
    
    C21 = C21(bc_data.unkn_d, bc_data.unkn_u);
    C22 = C22(bc_data.unkn_d, bc_data.unkn_v);
    
    M12u = M12(bc_data.unkn_d, bc_data.unkn_u);
    M12v = M12(bc_data.unkn_d, bc_data.unkn_v);
    
    %% Empty matrices
    uv0 = sparse(size(Au,1), size(Av,1));
    pp0 = sparse(size(G1,1), size(G1,1));
    pd0 = sparse(size(G1,1), size(B,1));
    
    %% Jacobian
    J = [                  Au                   uv0    T1  dt*theta*K21*Nr
                         uv0'                    Av    T2  dt*theta*K22*Nr
                           G1                    G2   pp0              pd0
         dt*theta*C21-M12u*Su  dt*theta*C22-M12v*Sv  pd0'               B];
end