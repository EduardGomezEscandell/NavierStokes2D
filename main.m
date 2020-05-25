clear variables;
clc; addpath('Post-processing');
%% Data entry

% Phisical data
width = 2;
height = 3;
duration = 5;

visc0 = 1;
mu  = 0.1;
omega = 1;

% Discretization
x_elems = 10;
y_elems = 15;

mesh.steps = 20;
theta = 1/2;

% Numerical parameters
method = 1;
% 1: Picard
% 2: Newton-Raphson

boundary_method = 1; 
% 1: row elimination
% 2: penalty method
penalty = 1e8;

stabilize_pressure = false;
stabilize_concentration = true;

% Precision
maxIter = 20;
tol = 1e-8;
relaxation = 0.3;

%% Generating mesh
Gamma_test{1} = @(X)(X(1) == 0 && X(2) < height/2 && X(2) > 0);
Gamma_test{2} = @(X)(X(1) == 0 && X(2) >= height/2  && X(2) <= height);
Gamma_test{3} = @(X)(X(2) == height && X(1) >0 && X(1) < width);
Gamma_test{4} = @(X)(X(1) == width);
Gamma_test{5} = @(X)(X(2) == 0 && X(1) >0 && X(1) < width);
[coords,connect, corner_to_node, node_to_corner, Gamma] = square_mesh(width, height, x_elems, y_elems, Gamma_test);

refelem.L1 = set_reference_element(1,1); % Linear 1D element
refelem.L2 = set_reference_element(1,2); % Quadratic 1D element
refelem.Q1 = set_reference_element(2,1); % Linear 2D quad element
refelem.Q2 = set_reference_element(2,2); % Quadraic 2D quad element

%% Shorthand notation
mesh.nodes = size(coords,2);
mesh.corners = (x_elems+1) * (y_elems+1);
mesh.elems = size(connect,2);
mesh.rows = y_elems;
mesh.cols = x_elems;

dof.u =               1:mesh.nodes;
dof.v = dof.u(end) + (1:mesh.nodes);
dof.p = dof.v(end) + (1:mesh.corners);
dof.d = dof.p(end) + (1:mesh.corners);

mesh.dof = max([dof.u, dof.v, dof.p, dof.d]);

dt = duration / mesh.steps;
%% Solution vectors
X_history = zeros(mesh.dof, mesh.steps+1);
Pe_history = zeros(mesh.nodes, mesh.steps+1);

%% Selecting boundary conditions
bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, dof, duration, omega);

%% Solver
% Constant terms
[K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, refelem, Gamma);

% Terms evaluated at step 1
X = X_history(:,1);
visc = get_viscosity(X(dof.p), visc0);
source = get_source(X(dof.u), X(dof.v));
[K1, ~, ~, C1, ~, ~, M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);
Pe_history(:,1) = get_peclet(coords, connect, mesh, dof, X, refelem, mu);

%% TEST STEADY STATE
Z1  = sparse(mesh.corners,mesh.corners);
Z12 = sparse(mesh.nodes,mesh.corners);
Z2  = sparse(mesh.nodes,mesh.nodes);

for i=1:maxIter
    visc = get_viscosity(X(dof.p), visc0);
    source = get_source(X(dof.u), X(dof.v));
    [K1, ~, ~, C1, ~, ~, M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);
    
%     h = - G1 * bc_data.U_dirichlet(:,1) + G2 * bc_data.V_dirichlet(:,1);
%     f1 = - M2*bc_data.U_dirichlet(:,1);
%     f2 = - M2*bc_data.V_dirichlet(:,1);
    
    B = -mu*K0 + C1;

%     B = B - C1*(K_tau + C1_tau);
    
    K = [ K1   Z2  G1' Z12
          Z2   K1  G2' Z12
          G1   G2  Z1   Z1
         Z12' Z12' Z1    B];

%     F = [f1; f2; h; M12*source];
    F = zeros(mesh.dof,1);
    
    F(dof.d) = (M12 - C1 * M12_tau)*source;
%     Boundary conditions
    K(bc_data.removed_rows,:) = sparse(bc_data.n_removed,mesh.dof);
    K = K + bc_data.dirichlet_matrix;
    F(bc_data.removed_rows) = 0;
    
    F(dof.u) = bc_data.U_dirichlet(:,1);
    F(dof.v) = bc_data.V_dirichlet(:,1);
    F(dof.p) = bc_data.P_neumann(:,1);
    F(dof.d) = bc_data.D_dirichlet(:,1) + M12*source;
    %
    X_new = K\F;
    error = norm(X-X_new)
    X = X + 1 *(X_new-X);
    post_processing_single(coords, connect, mesh, dof, corner_to_node, X)
    if  error < tol
        break;
    end
end
%% END OF TEST

for step = 1:mesh.steps
    %% Main iteration loop
    % Storing all terms that are a function of X^{n-1}
    C1_old = C1;
    K1_old = K1;
    source_old = source;
    
    for iter=1:maxIter
        %% Assembling elemental matrices
        source = get_source(X(dof.u), X(dof.v));
        visc = get_viscosity(X(dof.p), visc0);
        
        [K1, K21, K22, C1, C21, C22, supg] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);
        
        %% Assembling global system
        A = M2 + dt*theta*K1;
        T1 = dt*theta*G1';
        T2 = dt*theta*G2';
        b1 = (M2 + dt*(1-theta)*K1_old) * X(dof.u) + dt*(1-theta)*G1'*X(dof.p);
        b2 = (M2 + dt*(1-theta)*K1_old) * X(dof.v) + dt*(1-theta)*G2'*X(dof.p);
        
        B = M1 + dt*theta*(C1 + mu*K0);
        d = (M1 - dt*(1-theta)*(C1_old + mu*K0))*X(dof.d) + ...
                              dt*(theta*M12*source + (1-theta)*M12*source_old);
        
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
        
        K = [  A   Z2  -T1  Z12;
              Z2    A  -T2  Z12;
              G1   G2  Z1   Z1;
             Z12' Z12' Z1   B];

        F = [b1; b2; z; d];
% 
%         K = [  A   Z2  T1;
%               Z2    A  T2;
%               G1   G2  Z1];
% 
%         F = [b1; b2; z];

        
        %% Boundary condition enforcement
        if boundary_method==1
            K(bc_data.removed_dof,:) = sparse(bc_data.n_removed,mesh.dof);
            K = K + bc_data.dirichlet_matrix;
            F(bc_data.removed_dof) = 0;
            F = F + [bc_data.U_dirichlet(:,step);
                     bc_data.V_dirichlet(:,step);
                     bc_data.P_neumann(:,step)  ;
                     bc_data.D_dirichlet(:,step)];
        else
            % Nope
        end
        
        %% Obtaining residual and error
        
        R = F - K*X;
        error = max(R);

%         post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
%         fprintf('Error = %g\n',error);

        if error < tol
            break;
        end
        
        %% Calculating next X
        
        X_new = K\F;
        
        X = X + relaxation * (X_new-X);
        
    end
    %% Preparing next step
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    if(mod(step,10) == 0)
        fprintf('Step %3d of %3d completed\n',step, mesh.steps);
    end
    
    Pe_history(:,step+1) = get_peclet(coords, connect, mesh, dof, X, refelem, mu);
    X_history(:,step+1) = X;
end
fprintf('Step %3d of %3d completed\n',mesh.steps, mesh.steps);

post_processing(coords, X_history, Pe_history, duration, dof, mesh, corner_to_node)