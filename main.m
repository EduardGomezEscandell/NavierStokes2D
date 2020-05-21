clear variables;
clc; addpath('Post-processing');
%% Data entry

% Phisical data
width = 2;
height = 3;
duration = 5;

visc0 = 1;
mu = 1;
omega = 4;

% Discretization
x_elems = 10;
y_elems = 15;

mesh.steps = 100;
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
stabilize_concentration = false;

% Precision
maxIter = 10;
tol = 1e-8;
relaxation = 1;

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
mean_area = width*height / mesh.elems;
h = sqrt(mean_area);

%% Solution vectors
X_history = zeros(mesh.dof, mesh.steps+1);

%% Selecting boundary conditions
bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, dof, duration, omega);

%% Solver
% Constant terms
[K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, refelem, Gamma);

% Terms evaluated at step 1
X = X_history(1:dof.d-1,1);
visc = get_viscosity(X(dof.p), visc0);
source = get_source(X(dof.u), X(dof.v));
[K1,~,~, C1,~,~,~] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);


%% TEST STEADY STATE
% Z1  = sparse(mesh.corners,mesh.corners);
% Z12 = sparse(mesh.nodes,mesh.corners);
% Z2  = sparse(mesh.nodes,mesh.nodes);
% 
% K = [ K1  Z2 G1'
%       Z2  K1 G2'
%       G1 G2  Z1];
% 
% % f1 = M12' * bc_data.P_neumann(:,1);
% % f2 = M12' * bc_data.P_neumann(:,1);
% h = G1 * bc_data.U_dirichlet(:,1) + G2 * bc_data.V_dirichlet(:,1);
% % 
% F = [zeros(2*mesh.nodes,1); h];
% % F = zeros(mesh.dof - mesh.corners, 1);
% % 
% if boundary_method==1
%     K(bc_data.removed_dof,:) = sparse(bc_data.n_removed,mesh.dof-mesh.corners);
%     K = K + bc_data.dirichlet_matrix;
%     F(bc_data.removed_dof) = 0;
%     F = F + [bc_data.U_dirichlet(:,1); bc_data.V_dirichlet(:,1); bc_data.P_neumann(:,1)];
% else
%     bc_data.H = bc_data.H(:,1:(dof.d(1) - 1));
%     K = K + penalty * (bc_data.H'*bc_data.H);
%     F = F + penalty * bc_data.H'*(bc_data.e(:,1));
% end
% 
% X = K\F;

% post_processing_single(coords, connect, mesh, dof, corner_to_node, X) 
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
        
%         B = M1 + dt*theta*(C1 + mu*K0);
%         d = (M1 - dt*(1-theta)*(C1_old + mu*K0))*X(dof.d) + ...
%                               dt*(theta*M12*source + (1-theta)*M12*source_old);
        
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
        if boundary_method==1
            K(bc_data.removed_dof,:) = sparse(bc_data.n_removed,mesh.dof-mesh.corners);
            K = K + bc_data.dirichlet_matrix;
            F(bc_data.removed_dof) = 0;
            F = F + [bc_data.U_dirichlet(:,step);
                     bc_data.V_dirichlet(:,step);
                     bc_data.P_neumann(:,step)  ];
        else
            % Nope
        end
        
        %% Obtaining residual and error
        
        R = F - K*X;
        error = norm(R);

%         post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
%         fprintf('Error = %g\n',error);

%         [Pe_global, Pe] = get_peclet(coords, connect, X, linear_elem, mu);

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
    
    X_history(1:dof.d(1)-1,step+1) = X;
end
fprintf('Step %3d of %3d completed\n',mesh.steps, mesh.steps);

post_processing(coords, X_history, duration, dof, mesh, corner_to_node)