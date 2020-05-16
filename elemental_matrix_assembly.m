function [M, K, K1, K21, K22, C1, C21, C22, SUPG] = elemental_matrix_assembly(connect, coords, X, mesh, dof, visc, refelem, mu, theta, dt) 
        
        n_elems = size(connect,2);
        n_nodes = size(coords,2); 

        M   = sparse(mesh.nodes,mesh.nodes);
        SUPG= sparse(mesh.nodes,mesh.nodes);
        
        K   = sparse(mesh.nodes,mesh.nodes);
        K1  = sparse(mesh.nodes,mesh.nodes);
        K21 = sparse(mesh.nodes,mesh.nodes);
        K22 = sparse(mesh.nodes,mesh.nodes);
        
        C1  = sparse(mesh.nodes,mesh.nodes);
        C21 = sparse(mesh.nodes,mesh.nodes);
        C22 = sparse(mesh.nodes,mesh.nodes);
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            x = [X(dof.u(1) + nodes-1),  ...
                 X(dof.v(1) + nodes-1),  ...
                 X(dof.p(1) + nodes-1),  ...
                 X(dof.d(1) + nodes-1)];
                 v = visc(nodes);
                  
            local_mat = FEM_matrices(local_coordinates, refelem, x, v, mu, theta, dt);
            
            % Mass matrix
            M(nodes,nodes) 	 = M(nodes,nodes)   + local_mat.M;
            
            % Difusion matrices
            K(nodes,nodes) 	 = K(nodes,nodes)   + local_mat.K;
            K1(nodes,nodes)  = K1(nodes,nodes)  + local_mat.K1;
            K21(nodes,nodes) = K21(nodes,nodes) + local_mat.K21;
            K22(nodes,nodes) = K22(nodes,nodes) + local_mat.K22;
            
            % Convection matrices
            C1(nodes,nodes)  = C1(nodes,nodes)  + local_mat.C1;
            C21(nodes,nodes) = C21(nodes,nodes) + local_mat.C21;
            C22(nodes,nodes) = C22(nodes,nodes) + local_mat.C22;
            
            % Stabilization matrices
            SUPG(nodes,nodes) = SUPG(nodes,nodes) + local_mat.supg;
        end
end