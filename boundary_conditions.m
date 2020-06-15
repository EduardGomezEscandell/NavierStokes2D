function bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, duration)   
    %% Useful values
    span = max(coords');
    width = span(1);
    height= span(2);
    t = linspace(0, duration, mesh.steps+1);
    
    %% Initializing bc_data object
    bc_data = bc_data_init(mesh);
    
    %% Setting boundary conditions
    for node = [Gamma.nodes{1}, Gamma.nodes{2}]
        % Gamma 1 & 2
        bc_data = set_BC(bc_data, node, 'u', 0);
        bc_data = set_BC(bc_data, node, 'v', -1-sin(2*pi*t));
    end
    
    for node = Gamma.nodes{3}
        % Gamma 3
        bc_data = set_BC(bc_data, node, 'u', -1-sin(2*pi*t));
        bc_data = set_BC(bc_data, node, 'v', 0);
        bc_data = set_BC(bc_data, node, 'd', 1);
    end
    
    for node = Gamma.nodes{4}
        % Gamma 4
          bc_data = set_BC(bc_data, node, 'v', 1+sin(2*pi*t));
          bc_data = set_BC(bc_data, node, 'u', 0);
          bc_data = set_BC(bc_data, node, 'd', 0);
    end
    
    for node = Gamma.nodes{5}
        % Gamma 5
        bc_data = set_BC(bc_data, node, 'u', 1+sin(2*pi*t));
        bc_data = set_BC(bc_data, node, 'v', 0);
        bc_data = set_BC(bc_data, node, 'd', 0);
    end
    
    %% Writing pressure BC in case its enclosed flow
    
    if size(bc_data.enforced_p,1) == 0
        bc_data = set_BC(bc_data, 1, 'p', 0);
    end
    
    if size(bc_data.enforced_d,1) == 0
        bc_data = set_BC(bc_data, 1, 'd', 0);
    end
    
    %% Finishing bc_data object
    bc_data = bc_data_finish(bc_data, mesh);
    
    
    %% Assisting functions
    
    function bc_data = bc_data_init(mesh)
        bc_data.enforced_u = [];
        bc_data.enforced_v = [];
        bc_data.enforced_p = [];
        bc_data.enforced_d = [];
        
        bc_data.enforced_u_value = zeros(0, mesh.steps+1);
        bc_data.enforced_v_value = zeros(0, mesh.steps+1);
        bc_data.enforced_p_value = zeros(0, mesh.steps+1);
        bc_data.enforced_d_value = zeros(0, mesh.steps+1);
        
        bc_data.dirichlet_matrix = sparse(mesh.dof, mesh.dof);
    end

    function bc_data = bc_data_finish(bc_data, mesh)
        bc_data.unkn_u = setdiff(1:mesh.nodes, bc_data.enforced_u);
        bc_data.unkn_v = setdiff(1:mesh.nodes, bc_data.enforced_v);
        bc_data.unkn_p = setdiff(1:mesh.corners, bc_data.enforced_p);
        bc_data.unkn_d = setdiff(1:mesh.corners, bc_data.enforced_d);
    end
    
    function bc_data = set_BC(bc_data, node, variable, value)
        if(size(value,1) == 2)
           value = value * ones(1,mesh.steps)';
        end
        
        switch variable
            case 'u'
                bc_data.enforced_u_value(end+1,:) = value;
                bc_data.enforced_u(end+1) = node;
            case 'v'
                bc_data.enforced_v_value(end+1,:) = value;
                bc_data.enforced_v(end+1) = node;
            case 'p'
                corner = node_to_corner(node);
                if corner > 0
                    bc_data.enforced_p_value(end+1,:) = value;
                    bc_data.enforced_p(end+1) = corner;
                end
            case 'd'
                corner = node_to_corner(node);
                if corner > 0
                    bc_data.enforced_d_value(end+1,:) = value;
                    bc_data.enforced_d(end+1) = corner;
                end
            otherwise
                error('Unrecognized variable');
        end
    end
end





