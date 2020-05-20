function [K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, linear_elem, quadra_elem)
    % Linear
    M1 = sparse(mesh.corners,mesh.corners);
    
    % Quadratic
    K0  = sparse(mesh.corners,mesh.corners);
    M2 = sparse(mesh.nodes,mesh.nodes);
    
    % Hybrid
    M12 = sparse(mesh.corners,mesh.nodes);
    G1 = sparse(mesh.corners,mesh.nodes);
    G2 = sparse(mesh.corners,mesh.nodes);
    
    for el=1:mesh.elems
        nodes = connect(:,el);
        corners = node_to_corner(nodes(1:4));
        local_coords = coords(:,nodes);
        
        [m1, m12, m2, k, g1, g2] = FEM_constant(local_coords, linear_elem, quadra_elem);

        K0(corners,corners) = K0(corners,corners) + k;
        
        M1(corners,corners) = M1(corners,corners) + m1;
        M12(corners,nodes) =  M12(corners,nodes) + m12;
        M2(nodes,nodes)     = M2(nodes,nodes) + m2;
              
        G1(corners,nodes) = G1(corners,nodes) + g1;
        G2(corners,nodes) = G2(corners,nodes) + g2;
    end
end