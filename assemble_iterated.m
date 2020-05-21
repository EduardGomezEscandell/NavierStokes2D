function [K1, K21, K22, C1, C21, C22, Q, SUPG] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt) 

    SUPG= sparse(mesh.nodes,mesh.nodes);
    
    Q = sparse(mesh.nodes,mesh.nodes);
    
    K1  = sparse(mesh.nodes,mesh.nodes);
    K21 = sparse(mesh.nodes,mesh.nodes);
    K22 = sparse(mesh.nodes,mesh.nodes);

    C1  = sparse(mesh.corners,mesh.corners);
    C21 = sparse(mesh.corners,mesh.corners);
    C22 = sparse(mesh.corners,mesh.corners);

    for el=1:mesh.elems
        nodes = connect(:,el);
        corners = node_to_corner(nodes(1:4));
        local_coords = coords(:, nodes);

        vel = [X(dof.u(1) + nodes-1),  X(dof.v(1) + nodes-1)];
        conc = 0; %X(dof.d(1) + corners-1);
        v = visc(corners);

        local_mat = FEM_iterated(local_coords, vel, conc, v, mu, theta, dt, refelem.Q1, refelem.Q2);

        % Difusion matrices
        K1(nodes,nodes)  = K1(nodes,nodes)  + local_mat.K1;
        K21(nodes,nodes) = K21(nodes,nodes) + local_mat.K21;
        K22(nodes,nodes) = K22(nodes,nodes) + local_mat.K22;

        % Convection matrices
        C1(corners,corners)  = C1(corners,corners)  + local_mat.C1;
        C21(corners,corners) = C21(corners,corners) + local_mat.C21;
        C22(corners,corners) = C22(corners,corners) + local_mat.C22;

%             % Stabilization matrices
%             SUPG(corners,corners) = SUPG(corners,corners) + local_mat.supg;
    end

    for eg = 1:length(Gamma.neumann_edges)
        nodes = Gamma.neumann_edges(:,eg);
        corners = node_to_corner(nodes([1,3]));
        v = visc(corners);
        local_coords = coords(:,nodes);
        q = FEM_iterated_line(local_coords, v, refelem.L1, refelem.L2);
        Q(nodes, nodes) = Q(nodes, nodes) + q;
    end
end