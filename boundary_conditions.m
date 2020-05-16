function [free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, mesh, node_to_corner, dof, duration, omega)
    
    width=max(coords(1,:));
    height=max(coords(2,:));
    
    H = sparse(0,mesh.dof);
    e = sparse(0,mesh.steps+1);
    
    removed_dof = [];
    
    t = linspace(0, duration, mesh.steps+1);
    zero = zeros(size(t));
    one  = ones(size(t));
    sine = 1 + sin(omega*t - pi/2);
        
    for n = 1:mesh.nodes                   % Quadratic BC
        X = coords(:,n);
        
        dof_u = dof.u(1) + n - 1;
        dof_v = dof.v(1) + n - 1;
        dof_p = dof.p(1) + node_to_corner(n) - 1;
        dof_d = dof.d(1) + n - 1;
        
        if X(2) == 0 || X(2)==height
            % Gamma 5 and 3 (h walls)
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_v, zero);
            if node_to_corner(n) > 0 %( if n is a corner)
%                 [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_p, one);
            end
        elseif X(1) == 0 && X(2) < height/2
            % Gamma 1 (bottom left)
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_v, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_d, one);
            if node_to_corner(n) > 0 %( if n is a corner)
%                 [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_p, one);
            end
        elseif X(1) == 0 && X(2) >= height/2
            % Gamma 2 (top left)
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_v, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_d, one);
            if node_to_corner(n) > 0 %( if n is a corner)
%                 [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_p, one);
            end
        elseif X(1) == width
            % Gamma 4 (right wall)
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_u, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_v, zero);
            [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_d, zero);
            if node_to_corner(n) > 0 %( if n is a corner)
%                 [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof_p, one);
            end
        end
    end
    
    [removed_dof, X_history, H, e] = set_BC(removed_dof, X_history, H, e, dof.p(1), one);

    free_dof = 1:mesh.dof;
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


