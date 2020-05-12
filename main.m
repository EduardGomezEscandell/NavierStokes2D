clear variables;
clc;
%% Data entry

% Phisical data
width = 2;
height = 3;
duration = 5;

visc0 = 1;
mu = 1;
omega = 1;

% Discretization
x_elems = 20;
y_elems = 30;
degree = 1; %(space)

n_steps = 100;
theta = 1/2;

% Numerical parameters
method = 1;
% 1: Picard
% 2: Newton-Raphson

boundary_method = 1;
% 1: Penalty method
% 2: Row elimination
penalty = 1e8;

stabilize_pressure = false;
stabilize_concentration = false;

% Precision
maxIter = 50;
tol = 1e-5;
relaxation = 1;

%% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
refelem = set_reference_element(degree);

%% Shorthand notation
n_nodes = length(coords);
n_elems = length(connect);

u_dof = 1:n_nodes;
v_dof = n_nodes+1:2*n_nodes;
p_dof = 2*n_nodes+1:3*n_nodes;
d_dof = 3*n_nodes+1:4*n_nodes;

dt = duration / n_steps;
mean_area = width*height / n_elems;
h = sqrt(mean_area);

%% Solution vectors
X_history = zeros(n_nodes*4, n_steps);

%% Selecting boundary conditions
[free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, duration, omega);

%% Solver

for step = 2:n_steps
    %% Previous step info
    % Calculating all terms that are a function of X^{n-1}
    
    X = X_history(:,step-1);
    
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
    visc = visc0 + visc0 ./ (1 + exp(-10*(X(d_dof) - 0.5)));
        
    b1 = zeros(n_nodes,1);
    b2 = zeros(n_nodes,1);
    d_non_iter = zeros(n_nodes,1);
    Q1 = sparse(n_nodes,n_nodes);
    Q2 = sparse(n_nodes,n_nodes);
    
    for el=1:n_elems
        nodes = connect(:,el);
        local_coordinates = coords(:, nodes);
        
        x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*n_nodes), X(nodes+3*n_nodes)];
        v = visc(nodes);
        s = source(nodes);
        
        [M, K, K1, C1, q1, q2] = non_iterated_arrays(local_coordinates, x, v, refelem);
        
        b = (M + dt*(1-theta)*K1);
        d = (M - dt*(1-theta)*(C1 + mu*K));
        
        b1(nodes) = b1(nodes) + b * x(:,1) + dt * (1-theta)*q1*x(:,3);
        b2(nodes) = b2(nodes) + b * x(:,2) + dt * (1-theta)*q2*x(:,3);
        
        d_non_iter(nodes)  =  d_non_iter(nodes) ...
                  + d*x(:,4) + dt * (1-theta)*M*s;
              
        Q1(nodes,nodes) = Q1(nodes,nodes) + q1;
        Q2(nodes,nodes) = Q2(nodes,nodes) + q2;
    end
    
    %% Main iteration loop
    % Calculating all terms that are a function of X^{n}
    
    X = X_history(:,step);
    X(free_dof) = X_history(free_dof,step-1);
    
    for iter=1:maxIter
        %% Assembling elemental matrices
        modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
        source = 1./(1 + exp(-10*(modU - 0.5)));
        visc = visc0 + visc0 ./ (1 + exp(-10*(X(d_dof) - 0.5)));
        a = max(modU);
        
        [M, K, K1, K21, K22, C1, C21, C22, L, supg] = elemental_matrix_assembly(connect, coords, X, visc, refelem, mu, theta, dt);
        
        %% Assembling global system
        
        A = M + dt*theta*K1;
        T1 = -dt*theta*Q1;
        T2 = -dt*theta*Q2;
        
        B = M + dt*theta*(C1 + mu*K);
        d = d_non_iter + dt * theta * M * source;
        
        Z = zeros(size(A));
        z = zeros(size(d));
        
        if stabilize_concentration
            I = speye(size(B));
            B = (I + supg)*B;
            d = (I + supg)*d;
        end
        
        if ~stabilize_pressure
            L = Z;
        end
        
        K = [Q1  Q2   L  Z;
              A   Z  T1  Z;
              Z   A  T2  Z;
              Z   Z   Z  B];

        F = [ z; b1; b2; d];
        
        %% Boundary condition enforcement
        
        if boundary_method == 1             % Penalty method
            K = K + penalty * (H'*H);
            F = F + penalty * H'*e(:,step);
            system_size = 1:4*n_nodes;
        else                                % Row elimination
            K = K(free_dof,free_dof);
            F = F(free_dof);
            system_size = free_dof;
        end
        
        %% Obtaining residual and error
        
        R = K*X(system_size) - F;
        error = max(abs(R));

%         post_processing_single(coords, connect, X);
%         fprintf('Error = %g\n',error);

        [Pe_global, Pe] = get_peclet(coords, connect, X, refelem, mu);

        if error < tol
            break;
        end
        
        %% Calculating next X
        
        if method==1    % Picard
            X_new = K\F;
        else            % Newton-Rhapson
            % Calculating jacobian
            J = calc_jacobian(coords, visc0, theta, dt, X, K21, K22, C21, C22, Q1, Q2, T1, T2, M, A, B);
            J = J(system_size, system_size);

            % Newton Raphson
            X_new = X(system_size) - (J \ R);
        end
        
        X(system_size) = (1-relaxation)*X(system_size) + relaxation*X_new;
        
    end
    %% Preparing next step
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    if(mod(step,10) == 0)
        fprintf('Step %3d of %3d completed\n',step, n_steps);
    end
    
    X_history(:,step) = X;
end
fprintf('Step %3d of %3d completed\n',n_steps, n_steps);

post_processing(coords, X_history, duration)