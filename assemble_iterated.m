function [K1, K21, K22, C1, C21, C22, SUPG] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, visc, linear_elem, quadra_elem, mu, theta, dt) 

        SUPG= sparse(mesh.nodes,mesh.nodes);
        
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
                  
            local_mat = FEM_iterated(local_coords, vel, conc, v, mu, theta, dt, linear_elem, quadra_elem);
            
            % Difusion matrices
            K1(nodes,nodes)  = K1(nodes,nodes)  + local_mat.K1;
            K21(nodes,nodes) = K21(nodes,nodes) + local_mat.K21;
            K22(nodes,nodes) = K22(nodes,nodes) + local_mat.K22;
            
            % Convection matrices
            C1(corners,corners)  = C1(corners,corners)  + local_mat.C1;
            C21(corners,corners) = C21(corners,corners) + local_mat.C21;
            C22(corners,corners) = C22(corners,corners) + local_mat.C22;
            
            % Stabilization matrices
            SUPG(corners,corners) = SUPG(corners,corners) + local_mat.supg;
        end
end