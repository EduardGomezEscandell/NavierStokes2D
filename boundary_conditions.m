function [free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, duration, omega)
    n_steps = size(X_history, 2);
    n_nodes = size(coords, 2);
    width=max(coords(1,:));
    height=max(coords(2,:));
    
    H = sparse(0,4*n_nodes);
    e = sparse(0,n_steps);
    
    removed_dof = [];
    t = linspace(0, duration, n_steps);
    for n = 1:n_nodes
        X = coords(:,n);
        n_u = n;
        n_v = n + n_nodes;
        n_p = n + 2*n_nodes;
        n_r = n + 3*n_nodes;
        if X(2) == 0 || X(2)==height
            %% Gamma 3 and 5
            % U
            removed_dof(end+1) = n_u;
            X_history(n_u,:) = zeros(size(t));
            H(end+1,n_u) = 1;
            e(end+1,:) = zeros(size(t));

            % V
            removed_dof(end+1) = n_v;
            X_history(n_v,:) = zeros(size(t));
            H(end+1,n_v) = 1;
            e(end+1,:) = zeros(size(t));
        elseif X(1) == 0
            %% Gamma 1 and 2
            % U
            removed_dof(end+1) = n_u;
            X_history(n_u,:) = zeros(size(t));
            H(end+1,n_u) = 1;
            e(end+1,:) = zeros(size(t));

            % V
            removed_dof(end+1) = n_v;
            X_history(n_v,:) = 1 + sin(omega*t - pi/2);
            H(end+1,n_v) = 1;
            e(end+1,:) = 1 + sin(omega*t - pi/2);
            
            % P
            removed_dof(end+1) = n_p;
            X_history(n_p,:) = zeros(size(t));
            H(end+1,n_p) = 1;
            e(end+1,:) = zeros(size(t));
            
            % RHO
            removed_dof(end+1) = n_r;
            X_history(n_r,:) = ones(size(t));
            H(end+1,n_r) = 1;
            e(end+1,:) = ones(size(t));
        elseif X(1) == width
            %% Gamma 4
            % U
            removed_dof(end+1) = n_u;
            X_history(n_u,:) = zeros(size(t));
            H(end+1,n_u) = 1;
            e(end+1,:) = zeros(size(t));

            % V
            removed_dof(end+1) = n_v;
            X_history(n_v,:) = zeros(size(t));
            H(end+1,n_v) = 1;
            e(end+1,:) = zeros(size(t));
            
            % RHO
            removed_dof(end+1) = n_r;
            X_history(n_r,:) = zeros(size(t));
            H(end+1,n_r) = 1;
            e(end+1,:) = zeros(size(t));
        end
    end

    free_dof = 1:4*n_nodes;
    free_dof(removed_dof) = [];
end