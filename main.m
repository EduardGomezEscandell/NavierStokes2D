clear variables;
clc; addpath('Post-processing');
%% Data entry

% Phisical data
width = 2;
height = 3;
duration = 5;

visc0 = 100;
mu = 1;
omega = 1;

% Discretization
x_elems = 10;
y_elems = 15;

mesh.steps = 10;
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
[coords,connect, corner_to_node, node_to_corner] = square_mesh(width, height, x_elems, y_elems, 2);
linear_elem = set_reference_element(1);
quadra_elem = set_reference_element(2);

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
[removed_dof, dirichlet_matrix, H, e] =  boundary_conditions(coords, mesh, node_to_corner, dof, duration, omega);

%% Solver
% Constant terms
[K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, linear_elem, quadra_elem);

% Terms evaluated at step 1
X = X_history(1:dof.d-1,1);
visc = get_viscosity(X(dof.p), visc0);
source = get_source(X(dof.u), X(dof.v));
[K1,~,~, C1,~,~,~] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, visc, linear_elem, quadra_elem, mu, theta, dt);


%% TEST STEADY STATE
Z1  = sparse(mesh.corners,mesh.corners);
Z12 = sparse(mesh.nodes,mesh.corners);
Z2  = sparse(mesh.nodes,mesh.nodes);

u_dirichlet_id = removed_dof(removed_dof >= dof.u(1));
u_dirichlet_id = u_dirichlet_id(u_dirichlet_id <= dof.u(end));

v_dirichlet_id = removed_dof(removed_dof >= dof.v(1));
v_dirichlet_id = v_dirichlet_id(v_dirichlet_id <= dof.v(end));

U_dirichlet = zeros(mesh.nodes,1);
V_dirichlet = zeros(mesh.nodes,1);

U_dirichlet(u_dirichlet_id) = e(1:length(u_dirichlet_id),1);
V_dirichlet(v_dirichlet_id - mesh.nodes) = e(length(u_dirichlet_id)+(1:length(v_dirichlet_id)),1);

h = G1*U_dirichlet + G2*V_dirichlet;

K = [ K1  Z2  G1'
      Z2  K1  G2'
      G1  G2  Z1];
  
F = [zeros(2*mesh.nodes,1);
           h              ];




if boundary_method==1
    K(removed_dof,:) = sparse(length(removed_dof),mesh.dof-mesh.corners);
    K = K + dirichlet_matrix;
    F(removed_dof) = e(:,1);
else
    H = H(:,1:(dof.d(1) - 1));
    K = K + penalty * (H'*H);
    F = F + penalty*H'*(e(:,1));
end

X = K\F;

post_processing_single(coords, connect, mesh, dof, corner_to_node, X) 
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
        
        [K1, K21, K22, C1, C21, C22, supg] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, visc, linear_elem, quadra_elem, mu, theta, dt);
        
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
            K(removed_dof,:) = sparse(length(removed_dof),mesh.dof-mesh.corners);
            K = K + dirichlet_matrix;
            F(removed_dof) = e(:,step);
        else
            H = H(:,1:(dof.d(1) - 1));

            K = K + penalty * (H'*H);
            F = F - penalty*H'*(e(:,step));
        end
        
        %% Obtaining residual and error
        
        R = F - K*X;
        error = norm(R);

        post_processing_single(coords, connect, mesh, dof, corner_to_node, X);
        fprintf('Error = %g\n',error);

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
    
    X_history(:,step+1) = X;
end
fprintf('Step %3d of %3d completed\n',mesh.steps, mesh.steps);

post_processing(coords, X_history, duration, dof, mesh, corner_to_node)