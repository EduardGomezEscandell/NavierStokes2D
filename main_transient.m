clear variables;
clc; addpath('Post-processing');
%% Adjustable parameters

% Phisical data
width = 1;
height = 1;
duration = 2;

visc0 = 0.1;
mu   =  10;
omega = 1;

% Discretization
x_elems = 10;
y_elems = 10;

mesh.steps = 20;
theta = 2;

% Numerical parameters
solver = 2;
% 1: Picard
% 2: Newton-Raphson

% Relaxation
relaxation = 0.5; % Recommended 1 for NR and ~0.5 for Picard

% Precision
maxIter = 20;
tol = 1e-3;

%% Generating mesh
Gamma_test{1} = @(X)(X(1) == 0        && X(2) < height/2    && X(2) > 0);
Gamma_test{2} = @(X)(X(1) == 0        && X(2) >= height/2   && X(2) < height);
Gamma_test{3} = @(X)(X(2) == height);
Gamma_test{4} = @(X)(X(1) == width    && X(2) < height      && X(2) > 0);
Gamma_test{5} = @(X)(X(2) == 0);
[coords,connect, corner_to_node, node_to_corner, Gamma] = square_mesh(width, height, x_elems, y_elems, Gamma_test);

refelem.Q1 = set_reference_element(1); % Linear 2D quad element
refelem.Q2 = set_reference_element(2); % Quadraic 2D quad element

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

%% Selecting boundary conditions
bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, dof, duration, omega);

%% Solver
% Constant terms
[K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, refelem);

% Terms evaluated at step 1
X = X_history(:,1);
visc = get_viscosity(X(dof.p), visc0, 0);
source = get_source(X(dof.u), X(dof.v), 0);
[K1, ~, ~, C1, ~, ~, M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);
Pe_history(:,1) = get_peclet(coords, connect, mesh, dof, X, refelem, mu);

steady_state;

% Reducing matrices
G1red  = G1(bc_data.unkn_p, bc_data.unkn_u);
G2red  = G2(bc_data.unkn_p, bc_data.unkn_v);

for step = 1:mesh.steps
    %% Main iteration loop
    % Storing all terms that are a function of X^{n-1}
    C1_old = C1;
    K1_old = K1;
    source_old = source;
    
    % Setting up boundary conditions
    X(dof.u(1) -1 + bc_data.enforced_u) = bc_data.enforced_u_value(:,step+1);
    X(dof.v(1) -1 + bc_data.enforced_v) = bc_data.enforced_v_value(:,step+1);
    X(dof.p(1) -1 + bc_data.enforced_p) = bc_data.enforced_p_value(:,step+1);
    X(dof.d(1) -1 + bc_data.enforced_d) = bc_data.enforced_d_value(:,step+1);
    
    for iter=1:maxIter
        %% Assembling elemental matrices
        source = get_source(X(dof.u), X(dof.v), 0);
        visc = get_viscosity(X(dof.p), visc0, 0);
        
        [K1, K21, K22, C1, C21, C22,  M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);

        
        %% Assembling global system
        A = M2 + dt*theta*K1;
        T1 = dt*theta*G1';
        T2 = dt*theta*G2';
        b1 = (M2 + dt*(1-theta)*K1_old) * X(dof.u) + dt*(1-theta)*G1'*X(dof.p);
        b2 = (M2 + dt*(1-theta)*K1_old) * X(dof.v) + dt*(1-theta)*G2'*X(dof.p);
        
        B = M1 + dt*theta*(C1 + mu*K0);
        d = (M1 - dt*(1-theta)*(C1_old + mu*K0))*X(dof.d) + dt*(theta*M12*source + (1-theta)*M12*source_old);
        
        %% Reducing matrices
        Au = A(bc_data.unkn_u, bc_data.unkn_u);
        Av = A(bc_data.unkn_v, bc_data.unkn_v);
        T1red = T1(bc_data.unkn_u, bc_data.unkn_p);
        T2red = T2(bc_data.unkn_v, bc_data.unkn_p);
        Bred = B(bc_data.unkn_d, bc_data.unkn_d);
        
        %% Load vectors
        bu = - A(bc_data.unkn_u, bc_data.enforced_u)*bc_data.enforced_u_value(:,step+1) - T1(bc_data.unkn_u, bc_data.enforced_p)*bc_data.enforced_p_value(:,step+1);
        bv = - A(bc_data.unkn_v, bc_data.enforced_v)*bc_data.enforced_v_value(:,step+1) - T2(bc_data.unkn_v, bc_data.enforced_p)*bc_data.enforced_p_value(:,step+1);
        bp = - G1(bc_data.unkn_p, bc_data.enforced_u)*bc_data.enforced_u_value(:,step+1) - G2(bc_data.unkn_p, bc_data.enforced_v)*bc_data.enforced_v_value(:,step+1);
        dred = d(bc_data.unkn_d) - B(bc_data.unkn_d, bc_data.enforced_d) * bc_data.enforced_d_value(:,step+1);
        
        %% Obtaining residual
        [K, F] = build_monolithic(Au, Av, T1red, T2red, G1red, G2red, Bred, bu, bv, bp, dred);
        
        Xred = reduce_vec(X, bc_data, dof);
        res = F - K*Xred;

        %% Solver iteration
        if solver == 1
            % Picard
            dXred = K \ res;
        else
            % Newton-Raphson
            Jred = build_jacobian(mesh, dof, bc_data, X, dt, theta, visc0, Au, Av, Bred, M12, K21, K22, C21, C22, G1red, G2red, T1red, T2red);
            dXred = Jred \ res;
        end
        
        dX = enlarge_dX(dXred, bc_data, dof, mesh);
        
        %% Calculating next X        
        X = X + relaxation * dX;
        
        %% Evaluating error
        [error, ierror] = max(dX([dof.u, dof.v, dof.d])); % Ignoring pressure
        
%         post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
%         fprintf('Max error = %g at dof %g (%s)\n',error, ierror, errtype);
        
        if error < tol
            break;
        end
    end
    %% Preparing next step
    
    if(ierror < dof.v(1))
        errtype = 'velocity u';
    elseif(ierror < dof.p(1))
        errtype = 'velocity v';
    else
        errtype = 'concentration';
    end
    fprintf('Step %3d | %2d iters | Error %.3e | dof %5d (%s)\n',step, iter, error, ierror, errtype);
    
    X_history(:,step+1) = X;
end

post_processing(coords, X_history, Pe_history, duration, dof, mesh, corner_to_node)



function Xred = reduce_vec(X, bc_data, dof)
    u_dof = bc_data.unkn_u;
    v_dof = dof.v(1) - 1 + bc_data.unkn_v;
    p_dof = dof.p(1) - 1 + bc_data.unkn_p;
    d_dof = dof.d(1) - 1 + bc_data.unkn_d;
    
    Xred = X([u_dof, v_dof, p_dof, d_dof]);
end

function dX = enlarge_dX(dXred, bc_data, dof, mesh)
    u_dof = 1:length(bc_data.unkn_u);
    v_dof = u_dof(end) + (1:length(bc_data.unkn_v));
    p_dof = v_dof(end) + (1:length(bc_data.unkn_p));
    d_dof = p_dof(end) + (1:length(bc_data.unkn_d));
    
    dX = zeros(mesh.dof,1);
    
    dX(dof.u(1) -1 + bc_data.unkn_u) = dXred(u_dof);
    dX(dof.v(1) -1 + bc_data.unkn_v) = dXred(v_dof);
    dX(dof.p(1) -1 + bc_data.unkn_p) = dXred(p_dof);
    dX(dof.d(1) -1 + bc_data.unkn_d) = dXred(d_dof);
    
end


