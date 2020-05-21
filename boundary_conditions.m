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
        bc_data = set_BC(bc_data, node, 'p', sine);
%         bc_data = set_BC(bc_data, node, 'v', zero);
    end
    
    for node = Gamma.nodes{4}
%         bc_data = set_BC(bc_data, node, 'u', 5*ones);
%         bc_data = set_BC(bc_data, node, 'v', zero);
        bc_data = set_BC(bc_data, node, 'p', zeros);
    end
    
%     bc_data = set_BC(bc_data, 1, 'p', ones);
    
    I = speye(mesh.dof-mesh.corners,mesh.dof-mesh.corners);
    D = I;
    D(bc_data.removed_dof, bc_data.removed_dof) = 0;
    bc_data.dirichlet_matrix = I - D;
    bc_data.n_removed = length(bc_data.removed_dof);
    
    
    % Assisting functions
    
    
    function bc_data = bc_data_init(mesh)
        bc_data.U_dirichlet = sparse(mesh.nodes, mesh.steps);
        bc_data.V_dirichlet = sparse(mesh.nodes, mesh.steps);
        bc_data.P_neumann = sparse(mesh.corners, mesh.steps);
        bc_data.removed_dof = [];
    end
    
    function bc_data = set_BC(bc_data, node, variable, value)        
        dof_id = -1;
        switch variable
            case 'u'
                dof_id = dof.u(1) + node - 1;
                bc_data.U_dirichlet(node,:) = value;
            case 'v'
                dof_id = dof.v(1) + node - 1;
                bc_data.V_dirichlet(node,:) = value;
            case 'p'
                corner = node_to_corner(node);
                if corner > 0
                    dof_id = dof.p(1) + corner - 1;
                    bc_data.P_neumann(corner,:) = value;
                end
            case 'd'
                corner = node_to_corner(node);
                if corner > 0
                    dof_id = dof.d(1) + corner - 1;
                    bc_data.D_dirichlet(corner,:) = value;
                end
            otherwise
                error('Unrecognized variable');
        end
        
        % Row elimination
        if dof_id > 0
            bc_data.removed_dof(end+1) = dof_id;
        end

    %     % Penalty method
    %     H(end+1,dof) = 1;
    %     e(end+1,:) = value;
    end
end





