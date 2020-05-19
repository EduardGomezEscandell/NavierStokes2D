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
dof.d = dof.p(end) + (1:mesh.corners);

mesh.dof = max([dof.u, dof.v, dof.p, dof.d]);

dt = duration / mesh.steps;
mean_area = width*height / mesh.elems;
h = sqrt(mean_area);

%% Solution vectors
X_history = zeros(mesh.dof, mesh.steps+1);

%% Selecting boundary conditions
[X_history, H, e] =  boundary_conditions(coords, X_history, mesh, node_to_corner, dof, duration, omega);

%% Solver

for step = 2:mesh.steps+1
    %% Previous step info
    % Calculating all terms that are a function of X^{n-1}
    
    X = X_history(:,step-1);
    
    source = get_source(X(dof.u), X(dof.v));
    visc = get_viscosity(X(dof.d), visc0);
    
    b1 = zeros(mesh.nodes,1);
    b2 = zeros(mesh.nodes,1);
    d_non_iter = zeros(mesh.corners,1);
    
    K0  = sparse(mesh.corners,mesh.corners);
    M1 = sparse(mesh.corners,mesh.corners);
    M12 = sparse(mesh.corners,mesh.nodes);
    M2 = sparse(mesh.nodes,mesh.nodes);
    G1 = sparse(mesh.corners,mesh.nodes);
    G2 = sparse(mesh.corners,mesh.nodes);
    
    for el=1:mesh.elems
        nodes = connect(:,el);
        corners = node_to_corner(nodes(1:4));
        local_coords = coords(:,nodes);
        
        vel = [X(dof.u(1) + nodes-1),  X(dof.v(1) + nodes-1)];
        pres = X(dof.p(1) + corners-1);
        conc = X(dof.d(1) + corners-1);
        
        v = visc(corners);
        s = source(nodes);
        
        [m1, m12, m2, k, K1, C1, g1, g2] = non_iterated_arrays(local_coords, vel, v, linear_elem, quadra_elem);
        
        b = (m2 + dt*(1-theta)*K1);
        d = (m1 - dt*(1-theta)*(C1 + mu*k));
        
        b1(nodes) = b1(nodes) + b * vel(:,1) + dt * (1-theta)*g1'*pres;
        b2(nodes) = b2(nodes) + b * vel(:,2) + dt * (1-theta)*g2'*pres;
        
        d_non_iter(corners)  =  d_non_iter(corners) ...
                  + d*conc + dt * (1-theta)*m12*s;
        
        K0(corners,corners) = K0(corners,corners) + k;
        
        M1(corners,corners) = M1(corners,corners) + m1;
        M12(corners,nodes) =  M12(corners,nodes) + m12;
        M2(nodes,nodes)     = M2(nodes,nodes) + m2;
              
        G1(corners,nodes) = G1(corners,nodes) + g1;
        G2(corners,nodes) = G2(corners,nodes) + g2;
    end
    
    %% Main iteration loop
    % Calculating all terms that are a function of X^{n}
    
    X = X_history(:,step);
    
    for iter=1:maxIter
        %% Assembling elemental matrices
        source = get_source(X(dof.u), X(dof.v));
        visc = get_viscosity(X(dof.d), visc0);
        % a = max(modU);
        
        [K1, K21, K22, C1, C21, C22, supg] = elemental_matrix_assembly(connect, coords, node_to_corner, X, mesh, dof, visc, linear_elem, quadra_elem, mu, theta, dt);
        
        %% Assembling global system
        
        A = M2 + dt*theta*K1;
        T1 = dt*theta*G1';
        T2 = dt*theta*G2';
        
        B = M1 + dt*theta*(C1 + mu*K0);
        d = d_non_iter + dt * theta * M12 * source;
        
        % Empty matrices of appropiate sizes
        Z1  = sparse(mesh.corners,mesh.corners);
        Z12 = sparse(mesh.nodes,mesh.corners);
        Z2  = sparse(mesh.nodes,mesh.nodes);
        z = zeros(mesh.corners,1);
        
        if stabilize_concentration
            I = speye(size(B));
            B = (I + supg)*B;
            d = (I + supg)*d;
        end
        
        if ~stabilize_pressure
            %L = Z;
        end
        
%         K = [  A   Z2  T1  Z12;
%               Z2    A  T2  Z12;
%               G1   G2  Z1   Z1;
%              Z12' Z12' Z1   B];
% 
%         F = [b1; b2; z; d];

        K = [  A   Z2  T1;
              Z2    A  T2;
              G1   G2  Z1];

        F = [b1; b2; z];

        
        %% Boundary condition enforcement
        H = H(:,1:(dof.d(1) - 1));
        
        K = K + penalty * (H'*H);
        F = F - penalty*H'*(e(:,step));
        
        %% Obtaining residual and error
        
        R = F - K*X(1:(dof.d(1) - 1));
        error = norm(R);

        post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
        fprintf('Error = %g\n',error);

        [Pe_global, Pe] = get_peclet(coords, connect, X, linear_elem, mu);

        if error < tol
            break;
        end
        
        %% Calculating next X
        
        if method==1    % Picard
            X_new = K\F;
        else            % Newton-Rhapson
%             % Calculating jacobian
%             J = calc_jacobian(coords, visc0, theta, dt, X, K21, K22, C21, C22, G1, G2, T1, T2, M2, A, B);
% 
%             % Newton Raphson
%             X_new = X - (J \ R);
        end
        
        X = (1-relaxation)*X + relaxation*X_new;
        
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