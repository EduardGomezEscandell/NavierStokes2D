clear all;clc;

% Domain
width = 2;
height = 3;
duration = 1;

% Discretization
x_elems = 8;
y_elems = 12;
n_steps = 10;

degree = 1; %(space)

% Non-linearity
maxIter = 50;
tol = 1e-5;

% Physical properties
visc0 = 1;
mu = 1;
omega = 1;

% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
refelem = set_reference_element(degree);

n_nodes = length(coords);
n_elems = length(connect);
nodes_per_elem = size(connect,1);

% Selecting boundary conditions
removed_dof = [];
enforced_vals = [];
t = linspace(0, duration, n_steps);
for n = 1:n_nodes
    X = coords(:,n);
    if X(1) == 0 % Gamma 1 and 2
        removed_dof(end+1) = n;
        enforced_vals(end+1,:) = 1 + sin(omega*t - pi/2);
        
        removed_dof(end+1) = n + 2*n_nodes;
        enforced_vals(end+1,:) = zeros(size(t));
        
        removed_dof(end+1) = n + 2*n_nodes;
        enforced_vals(end+1,:) = 1 * ones(1,n_steps);
        
    elseif X(2) == 0 % Gamma 1 and 2
        removed_dof(end+1) = n + 2*n_nodes;
        enforced_vals(end+1,:) = 0 * ones(1,n_steps);
    end
end

free_dof = 1:3*n_nodes;
free_dof(removed_dof) = [];

% Solution vectors
X_history = zeros(n_nodes*3, n_steps);
X_history(removed_dof,:) = enforced_vals;

% Initializing
dt = duration / n_steps;

% Assembling

for step = 2:n_steps
    X = X_history(:,step-1);
    X(removed_dof) = enforced_vals(:,step) - enforced_vals(:,step-1);   
    dX = X - X_history(:,step-1);
    
    F_global_non_iter = sparse(n_nodes*3,1);
    
    modU = sqrt(X(1:n_nodes).^2 + X(n_nodes+1:n_nodes*2).^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
    rho = X(2*n_nodes+1:3*n_nodes);
    visc = visc0 + visc0 ./ (1 + exp(-10*(rho - 0.5)));
        
    % Obtaining arrays and matrices that need not iteration
    for el=1:n_elems
        nodes = connect(:,el);
        local_coordinates = coords(:, nodes);
        
        x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
        dx = [dX(nodes),dX(nodes+n_nodes),dX(nodes+2*nodes)];
        v = visc(nodes);
        
        [K, K1, C1, M] = non_iterated_arrays(local_coordinates, x, dx, v, refelem);
        
        z = zeros(size(K1));
        b(1:nodes_per_elem,1)   = K1 * x(:,1);
        b(nodes_per_elem+1:nodes_per_elem*2,1) = K1 * x(:,2);
        d = -(C1 + mu*K)*x(:,3) + 0.5 * M*source(nodes);
        
        F_local = [b; d];
        
        dof = [nodes; nodes+n_nodes; nodes+2*n_nodes];
        
        F_global_non_iter(dof) = F_global_non_iter(dof) + F_local;
    end
    
    % Iterating
    K_global = sparse(n_nodes*3, n_nodes*3); % ux, uy, p, r
    
    
    for iter=1:maxIter
        k_non_zero = nnz(K_global);
        K_global = spalloc(n_nodes*3, n_nodes*3, k_non_zero);
        
        F_global = F_global_non_iter;
        
        modU = sqrt(X(1:n_nodes).^2 + X(n_nodes+1:n_nodes*2).^2);
        source = 1./(1 + exp(-10*(modU - 0.5)));
        rho = X(2*n_nodes+1:3*n_nodes);
        visc = visc0 + visc0 ./ (1 + exp(-10*(rho - 0.5)));
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
            dx = [dX(nodes),dX(nodes+n_nodes),dX(nodes+2*nodes)];
            s = source(nodes);
            v = visc(nodes);
            
            [K, K1, K21, K22, C1, C21, C22, M] = FEM_matrices(local_coordinates,refelem, x, dx, v);
            
            % Assembly of Solvent
            % A*u = 0
            a = 1/dt*M + 0.5*K1;
            z = zeros(size(K1));
            A =  [a z
                  z a];

            % Assembly of solute
            % B*u = S
            B = (1/dt*M + 0.5*C1 + mu*K);
            d = 0.5 * M * s;
            
            % Joining two systems
            % K*u = F
            K_local = [                                     A, zeros(2*nodes_per_elem,nodes_per_elem);
                       zeros(nodes_per_elem,2*nodes_per_elem),                                      B];
                   
            F_local = [zeros(2*nodes_per_elem,1); d];
                      
            % Assembling
            dof = [nodes;
                  nodes + n_nodes;
                  nodes + n_nodes*2];

            K_global(dof,dof) = K_global(dof, dof) + K_local;
            F_global(dof) = F_global(dof) + F_local; 
        end
        
        % Solving reduced system
        K_global = K_global(free_dof, free_dof);
        F_global = F_global(free_dof);
        
        % Newton-Rhapson
        R = K_global*dX(free_dof) - F_global;
        
        % Calculating jacobian
        J = calc_jacobian(coords, connect, X, dX, visc0, mu, dt, removed_dof,refelem);
        J = J(free_dof, free_dof);
        
        % Newton Rhapson
        dX(free_dof) = dX(free_dof) - J \ R;
        
        % Convergence condition
        R = K_global*dX(free_dof) - F_global;
        error = norm(R);
        if error < tol
            break;
        end
         X = X_history(:,step-1) + dX;
    end
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    Error(step) = error;
    
    if(mod(step - 1,10) == 0)
        fprintf('Step %3d of %3d completed\n',step, n_steps);
    end
    
    X_history(:,step) = X_history(:,step-1) + dX;
end

post_porcessing(coords, X_history, duration)