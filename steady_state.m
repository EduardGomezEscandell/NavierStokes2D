G1red  = G1(bc_data.unkn_p, bc_data.unkn_u);
G2red  = G2(bc_data.unkn_p, bc_data.unkn_v);

X(dof.u(1) -1 + bc_data.enforced_u) = bc_data.enforced_u_value(:,1);
X(dof.v(1) -1 + bc_data.enforced_v) = bc_data.enforced_v_value(:,1);
X(dof.p(1) -1 + bc_data.enforced_p) = bc_data.enforced_p_value(:,1);
X(dof.d(1) -1 + bc_data.enforced_d) = bc_data.enforced_d_value(:,1);


for i=1:maxIter
    visc = get_viscosity(X(dof.p), visc0);
    source = get_source(X(dof.u), X(dof.v));
    [K1, ~, ~, C1, ~, ~, M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);

    B = -mu*K0 + C1;
    d = (M12 - C1 * M12_tau)*source;
    
    % Matrices eliminations
    K1u = K1(bc_data.unkn_u, bc_data.unkn_u);
    K1v = K1(bc_data.unkn_v, bc_data.unkn_v);    
    Bred = B(bc_data.unkn_d, bc_data.unkn_d);
    
    % Load vectors
    bu = - K1(bc_data.unkn_u, bc_data.enforced_u)*bc_data.enforced_u_value(:,1) - G1(bc_data.enforced_p, bc_data.unkn_u)'*bc_data.enforced_p_value(:,1);
    bv = - K1(bc_data.unkn_v, bc_data.enforced_v)*bc_data.enforced_v_value(:,1) - G2(bc_data.enforced_p, bc_data.unkn_v)'*bc_data.enforced_p_value(:,1);
    bp = - G1(bc_data.unkn_p, bc_data.enforced_u)*bc_data.enforced_u_value(:,1) - G2(bc_data.unkn_p, bc_data.enforced_v) *bc_data.enforced_v_value(:,1);
    dred = d(bc_data.unkn_d) - B(bc_data.unkn_d, bc_data.enforced_d) * bc_data.enforced_d_value(:,1);
    
    [K, F] = build_monolythic(K1u, K1v, G1red, G2red, Bred, dred, bu, bv, bp);
    
    Xred = reduce_vec(X, bc_data, dof);
    res = F - K*Xred;
    
    dXred = K\res;
    dX = build_dX(dXred, bc_data, dof, mesh);
    
    error = norm(dX([dof.u, dof.v, dof.d])); % Ignoring pressure
    X = X + dX;
    
    post_processing_single(coords, connect, mesh, dof, corner_to_node, X)
    
    if  error < tol
        break;
    end
end

post_processing_single(coords, connect, mesh, dof, corner_to_node, X)


function dX = build_dX(dXred, bc_data, dof, mesh)   
    u_dof = 1:length(bc_data.unkn_u);
    v_dof = u_dof(end) + (1:length(bc_data.unkn_v));
    p_dof = v_dof(end) + (1:length(bc_data.unkn_p));
    d_dof = p_dof(end) + (1:length(bc_data.unkn_d));
    
    dX = zeros(mesh.dof,1);
    
    dX(dof.u(1) -1 + bc_data.unkn_u) = dXred(u_dof);
    dX(dof.v(1) -1 + bc_data.unkn_v) = dXred(v_dof);
    dX(dof.p(1) -1 + bc_data.unkn_p) = dXred(p_dof);
    dX(dof.d(1) -1 + bc_data.unkn_d) = dXred(d_dof);
    
end

function Xred = reduce_vec(X, bc_data, dof)
    u_dof = bc_data.unkn_u;
    v_dof = dof.v(1) - 1 + bc_data.unkn_v;
    p_dof = dof.p(1) - 1 + bc_data.unkn_p;
    d_dof = dof.d(1) - 1 + bc_data.unkn_d;
    
    Xred = X([u_dof, v_dof, p_dof, d_dof]);
end

function [K, F] = build_monolythic(K1u, K1v, G1, G2, B, d, bu, bv, bp)
    zuv = sparse(size(K1u,1), size(K1v,1));
    zpp = sparse(size(G1,1), size(G1,1));
    zud = sparse(size(K1u,1), size(B,1));
    zvd = sparse(size(K1v,1), size(B,1));
    zpd = sparse(size(G1,1), size(B,1));
    
    K = [ K1u    zuv    G1'     zud
          zuv'   K1v    G2'     zvd
          G1     G2     zpp     zpd
          zud'   zvd'   zpd'      B  ];
      
      
    F = [bu;
         bv;
         bp;
          d];
end