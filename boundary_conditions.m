function bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, dof, duration, omega)
    
    width=max(coords(1,:));
    height=max(coords(2,:));

    t = linspace(0, duration, mesh.steps);
    
    % Different arrays to make notation cleaner
    zero = zeros(size(t));
    one  = ones(size(t));
    sine = 2*sin(omega*t - pi/2);
    parabh = @(X)((X(2) - height/2)*(height-X(2)) * 16 / height^2);
    parabf = @(X)(X(2)*(X(2)-height) * 4 / height^2);
    
    %
    bc_data = bc_data_init(mesh);
    
    for node = [Gamma.nodes{1}, Gamma.nodes{3}, Gamma.nodes{5}]
        bc_data = set_BC(bc_data, node, 'u', zero);
        bc_data = set_BC(bc_data, node, 'v', zero);
    end
    
    for node = Gamma.nodes{2}
%         bc_data = set_BC(bc_data, node, 'u', 10*ones);
        bc_data = set_BC(bc_data, node, 'u', -sine);
%         bc_data = set_BC(bc_data, node, 'v', zero);
        bc_data = set_BC(bc_data, node, 'd', 1*ones);
    end
    
    for node = Gamma.nodes{4}
%         bc_data = set_BC(bc_data, node, 'u', 5*ones);
%         bc_data = set_BC(bc_data, node, 'v', zero);
        bc_data = set_BC(bc_data, node, 'p', ones);
        bc_data = set_BC(bc_data, node, 'd', 0*ones);
    end
    
%     bc_data = set_BC(bc_data, 1, 'p', ones);
    
    bc_data.n_removed = length(bc_data.removed_rows);
    
    % Assisting functions
    
    
    function bc_data = bc_data_init(mesh)
        bc_data.U_dirichlet = sparse(mesh.nodes, mesh.steps);
        bc_data.V_dirichlet = sparse(mesh.nodes, mesh.steps);
        bc_data.P_neumann = sparse(mesh.corners, mesh.steps);
        bc_data.D_dirichlet = sparse(mesh.corners, mesh.corners);
        bc_data.removed_rows = [];
        bc_data.dirichlet_matrix = sparse(mesh.dof, mesh.dof);
    end
    
    function bc_data = set_BC(bc_data, node, variable, value)        
        dof_id = -1;
        switch variable
            case 'u'
                dof_id = dof.u(1) + node - 1;
                bc_data.U_dirichlet(node,:) = value;
                removed_row = dof_id;
                bc_data.dirichlet_matrix(removed_row,dof_id) = 1; 
                bc_data.removed_rows(end+1) = removed_row;
            case 'v'
                dof_id = dof.v(1) + node - 1;
                bc_data.V_dirichlet(node,:) = value;
                removed_row = dof_id;
                bc_data.dirichlet_matrix(removed_row,dof_id) = 1;
                bc_data.removed_rows(end+1) = removed_row; 
            case 'p'
                corner = node_to_corner(node);
                if corner > 0
                    dof_id = dof.p(1) + corner - 1;
                    bc_data.P_neumann(corner,:) = value;
                    removed_row = dof.p(1) + corner - 1;
                    bc_data.dirichlet_matrix(removed_row,dof_id) = 1; 
                    bc_data.removed_rows(end+1) = removed_row;
                end
            case 'd'
                corner = node_to_corner(node);
                if corner > 0
                    dof_id = dof.d(1) + corner - 1;
                    bc_data.D_dirichlet(corner,:) = value;
                    removed_row = dof.d(1) + corner - 1;
                    bc_data.dirichlet_matrix(removed_row,dof_id) = 1; 
                    bc_data.removed_rows(end+1) = removed_row;
                end
            otherwise
                error('Unrecognized variable');
        end
    end
end





