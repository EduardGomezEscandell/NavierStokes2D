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
        % U
        removed_dof(end+1) = n;
        enforced_vals(end+1,:) = 1 + sin(omega*t - pi/2);
        
        % V
        removed_dof(end+1) = n + n_nodes;
        enforced_vals(end+1,:) = zeros(size(t));
        
        % RHO
        removed_dof(end+1) = n + 2*n_nodes;
        enforced_vals(end+1,:) = 1 * ones(1,n_steps);
        
    elseif X(2) == 0 % Gamma 4
        % RHO
        removed_dof(end+1) = n + 2*n_nodes;
        enforced_vals(end+1,:) = 0 * ones(1,n_steps);
    end
end

free_dof = 1:3*n_nodes;
free_dof(removed_dof) = [];

% Shorthand notation
u_dof = 1:n_nodes;
v_dof = n_nodes+1:2*n_nodes;
p_dof = 2*n_nodes+1:3*n_nodes;

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
    
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
    visc = visc0 + visc0 ./ (1 + exp(-10*(X(p_dof) - 0.5)));
        
    % Obtaining arrays and matrices that need not iteration
    b1 = zeros(n_nodes,1);
    b2 = zeros(n_nodes,1);
    d_non_iter = zeros(n_nodes,1);
    
    for el=1:n_elems
        nodes = connect(:,el);
        local_coordinates = coords(:, nodes);
        
        x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
        dx = [dX(nodes),dX(nodes+n_nodes),dX(nodes+2*nodes)];
        v = visc(nodes);
        
        [K, K1, C1, M] = non_iterated_arrays(local_coordinates, x, dx, v, refelem);
        
        z = zeros(size(K1));
        b1(nodes) = b1(nodes) + K1 * x(:,1);
        b2(nodes) = b2(nodes) + K1 * x(:,2);
        d_non_iter(nodes)  =  d_non_iter(nodes) - (C1 + mu*K)*x(:,3) + 0.5 * M*source(nodes);
    end
    
    % Iterating
    for iter=1:maxIter
        
        K   = sparse(n_nodes,n_nodes);
        M   = sparse(n_nodes,n_nodes);
        K1  = sparse(n_nodes,n_nodes);
        K21 = sparse(n_nodes,n_nodes);
        K22 = sparse(n_nodes,n_nodes);
        C1  = sparse(n_nodes,n_nodes);
        C21 = sparse(n_nodes,n_nodes);
        C22 = sparse(n_nodes,n_nodes);
        
        modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
        source = 1./(1 + exp(-10*(modU - 0.5)));
        visc = visc0 + visc0 ./ (1 + exp(-10*(X(p_dof) - 0.5)));
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
            dx = [dX(nodes),dX(nodes+n_nodes),dX(nodes+2*nodes)];
            s = source(nodes);
            v = visc(nodes);
                  
            % Assembling
            local_mat = FEM_matrices(local_coordinates, refelem, x, dx, v);
        
            K(nodes,nodes) 	 = K(nodes,nodes)   + local_mat.K;
            M(nodes,nodes) 	 = M(nodes,nodes)   + local_mat.M;
            K1(nodes,nodes)  = K1(nodes,nodes)  + local_mat.K1;
            K21(nodes,nodes) = K21(nodes,nodes) + local_mat.K21;
            K22(nodes,nodes) = K22(nodes,nodes) + local_mat.K22;
            C1(nodes,nodes)  = C1(nodes,nodes)  + local_mat.C1;
            C21(nodes,nodes) = C21(nodes,nodes) + local_mat.C21;
            C22(nodes,nodes) = C22(nodes,nodes) + local_mat.C22;
        end
        
        A = 1/dt * M + 0.5 * K1;
        B = 1/dt * M + 0.5*(C1 + mu*K);
        d = d_non_iter + 0.5*M*source;
        
        % Newton-Rhapson
        
        R(u_dof,:) = A*X(u_dof) - b1;
        R(v_dof,:) = A*X(v_dof) - b1;
        R(p_dof,:) = B*X(p_dof) - d;
        
        R = R(free_dof);
        
        % Convergence condition
        error = norm(R)
        if error < tol
            break;
        end
        
        % Calculating jacobian
        J = calc_jacobian(coords, visc0, X, K21, K22, C21, C22, M, A, B);
        J = J(free_dof, free_dof);
        
        % Newton Raphson
        dX(free_dof) = dX(free_dof) - J \ R;
        X(free_dof) = X_history(free_dof,step-1) + dX(free_dof);
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

post_processing(coords, X_history, duration)