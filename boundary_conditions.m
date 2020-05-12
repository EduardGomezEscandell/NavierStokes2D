function [free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, duration, omega)
    
    n_steps = size(X_history, 2);
    n_nodes = size(coords, 2);
    
    width=max(coords(1,:));
    height=max(coords(2,:));
    
    H = sparse(0,4*n_nodes);
    e = sparse(0,n_steps);
    
    removed_dof = [];
    
    t = linspace(0, duration, n_steps);
    zero = zeros(size(t));
    one  = ones(size(t));
    sine = 1 + sin(omega*t - pi/2);
    
    for n = 1:n_nodes
        X = coords(:,n);
        
        % dof shorthands:
        n_u = n;
        n_v = n + n_nodes;
        n_p = n + 2*n_nodes;
        n_d = n + 3*n_nodes;
        
        if X(2) == 0 || X(2)==height
            % Gamma 3 and 5
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_v, zero);
        elseif X(1) == 0
            % Gamma 1 and 2
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_v, sine);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_p, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_d, one);
        elseif X(1) == width
            % Gamma 4
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_v, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, n_d, zero);
        end
    end

    free_dof = 1:4*n_nodes;
    free_dof(removed_dof) = [];
end


function [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof, value)
    X_history(dof,:) = value;

    % Row elimination
    removed_dof(end+1) = dof;
    
    % Penalty method
    H(end+1,dof) = 1;
    e(end+1,:) = value;
end


