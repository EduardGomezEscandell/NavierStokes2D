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
x_elems = 8;
y_elems = 12;

mesh.steps = 10;
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
maxIter = 10;
tol = 1e-8;
relaxation = 1;

%% Generating mesh
[coords,connect, corner_to_node, node_to_corner] = square_mesh(width, height, x_elems, y_elems, 2);
linear_elem = set_reference_element(1);
quadra_elem = set_reference_element(2);

%% Shorthand notation
mesh.nodes = size(coords,2);
mesh.corners = (x_elems+1) * (y_elems+1);
mesh.elems = size(connect,2);

dof.u = 1:mesh.nodes;
dof.v = dof.u(end) + (1:mesh.nodes);
dof.p = dof.v(end) + (1:mesh.corners);
dof.d = dof.p(end) + (1:mesh.nodes);

mesh.dof = max([dof.u, dof.v, dof.p, dof.d]);

dt = duration / mesh.steps;
mean_area = width*height / mesh.elems;
h = sqrt(mean_area);

%% Solution vectors
X_history = zeros(mesh.dof, mesh.steps+1);

%% Selecting boundary conditions
[free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, mesh, node_to_corner, dof, duration, omega);

%% Solver

for step = 2:mesh.steps+1
    %% Previous step info
    % Calculating all terms that are a function of X^{n-1}
    
    X = X_history(:,step-1);
    
    modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
    visc = visc0 + visc0 ./ (1 + exp(-10*(X(dof.d) - 0.5)));
        
    b1 = zeros(mesh.nodes,1);
    b2 = zeros(mesh.nodes,1);
    d_non_iter = zeros(mesh.nodes,1);
    G1 = sparse(mesh.corners,mesh.nodes);
    G2 = sparse(mesh.corners,mesh.nodes);
    
    for el=1:mesh.elems
        nodes = connect(:,el);
        corners = node_to_corner(nodes(1:4));
        local_coordinates = coords(:, nodes);
        
       x = [X(dof.u(1) + nodes-1),  ...
            X(dof.v(1) + nodes-1),  ...
            X(dof.p(1) + nodes-1),  ...
            X(dof.d(1) + nodes-1)]; 

        v = visc(nodes);
        s = source(nodes);
        
        [M, K, K1, C1, g1, g2] = non_iterated_arrays(local_coordinates, x, v, linear_elem, quadra_elem);
        
        b = (M + dt*(1-theta)*K1);
        d = (M - dt*(1-theta)*(C1 + mu*K));
        
        b1(nodes) = b1(nodes) + b * x(:,1) + dt * (1-theta)*g1'*x(1:4,3);
        b2(nodes) = b2(nodes) + b * x(:,2) + dt * (1-theta)*g2'*x(1:4,3);
        
        d_non_iter(nodes)  =  d_non_iter(nodes) ...
                  + d*x(:,4) + dt * (1-theta)*M*s;
        
        G1(corners,nodes) = G1(corners,nodes) + g1;
        G2(corners,nodes) = G2(corners,nodes) + g2;
    end
    
    %% Main iteration loop
    % Calculating all terms that are a function of X^{n}
    
    X = X_history(:,step);
    
    for iter=1:maxIter
        %% Assembling elemental matrices
        modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
        source = 1./(1 + exp(-10*(modU - 0.5)));
        visc = visc0 + visc0 ./ (1 + exp(-10*(X(dof.d) - 0.5)));
        a = max(modU);
        
        [M, K, K1, K21, K22, C1, C21, C22, supg] = elemental_matrix_assembly(connect, coords, X, mesh, dof, visc, quadra_elem, mu, theta, dt);
        
        %% Assembling global system
        
        A = M + dt*theta*K1;
        T1 = -dt*theta*G1';
        T2 = -dt*theta*G2';
        
        B = M + dt*theta*(C1 + mu*K);
        d = d_non_iter + dt * theta * M * source;
        
        Z1 = sparse(mesh.corners,mesh.corners);
        Z12 = sparse(mesh.nodes,mesh.corners);
        Z2 = sparse(mesh.nodes,mesh.nodes);
        z = zeros(mesh.corners,1);
        
        if stabilize_concentration
            I = speye(size(B));
            B = (I + supg)*B;
            d = (I + supg)*d;
        end
        
        if ~stabilize_pressure
            %L = Z;
        end
        
        K = [ A  Z2  T1  Z2;
             Z2   A  T2  Z2;
             G1  G2  Z1  Z12';
             Z2  Z2  Z12  B];

        F = [b1; b2; z; d];
        
        %% Boundary condition enforcement
        
        if boundary_method == 1             % Penalty method
            K = K + penalty * (H'*H);
            F = F + penalty * H'*e(:,step);
            system_size = 1:mesh.dof;
        else                                % Row elimination
            K = K(free_dof,free_dof);
            F = F(free_dof);
            system_size = free_dof;
        end
        
        %% Obtaining residual and error
        
        R = K*X(system_size) - F;
        error = norm(R);

%         post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
%         fprintf('Error = %g\n',error);

        [Pe_global, Pe] = get_peclet(coords, connect, X, linear_elem, mu);

        if error < tol
            break;
        end
        
        %% Calculating next X
        
        if method==1    % Picard
            X_new = K\F;
        else            % Newton-Rhapson
            % Calculating jacobian
            J = calc_jacobian(coords, visc0, theta, dt, X, K21, K22, C21, C22, G1, G2, T1, T2, M, A, B);
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
        fprintf('Step %3d of %3d completed\n',step, mesh.steps);
    end
    
    X_history(:,step) = X;
end
fprintf('Step %3d of %3d completed\n',mesh.steps, mesh.steps);

post_processing(coords, X_history, duration, dof, mesh, corner_to_node)