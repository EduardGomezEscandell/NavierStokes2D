function [M, K, K1, K21, K22, C1, C21, C22, L, SUPG] = elemental_matrix_assembly(connect, coords, X, visc, refelem, mu, theta, dt) 
        
        n_elems = size(connect,2);
        n_nodes = size(coords,2); 

        M   = sparse(n_nodes,n_nodes);
        SUPG= sparse(n_nodes,n_nodes);
        
        K   = sparse(n_nodes,n_nodes);
        L   = sparse(n_nodes,n_nodes);
        K1  = sparse(n_nodes,n_nodes);
        K21 = sparse(n_nodes,n_nodes);
        K22 = sparse(n_nodes,n_nodes);
        
        C1  = sparse(n_nodes,n_nodes);
        C21 = sparse(n_nodes,n_nodes);
        C22 = sparse(n_nodes,n_nodes);
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            x  = [ X(nodes), X(nodes+n_nodes), X(nodes+2*nodes), X(nodes+3*nodes)];
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
            L(nodes,nodes) = L(nodes,nodes) + local_mat.L;
            SUPG(nodes,nodes) = SUPG(nodes,nodes) + local_mat.supg;
        end
end