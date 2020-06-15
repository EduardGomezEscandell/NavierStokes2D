clear variables;
clc; addpath('Post-processing');
%% Adjustable parameters

% Physical data
width = 1;
height = 1;

visc0 = 0.1;
mu   =  10;
omega = 1;

% Discretization
x_elems = 10;
y_elems = 10;

% Numerical parameters
solver = 2;
% 1: Picard
% 2: Newton-Raphson

% Relaxation
relaxation = 0.1; % Recommended 1 for NR and ~0.5 for Picard

% Precision
maxIter = 20;
tol = 1e-8;

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

%% Dummy values
mesh.steps = 1;
dt = 0;
theta = 0;
duration = 0;

%% Selecting boundary conditions
bc_data =  boundary_conditions(coords, mesh, Gamma, node_to_corner, dof, duration, omega);

%% Solution vectors
X = zeros(mesh.dof, 1);

%% Solver
% Constant terms
[K0, M1, M12, M2, G1, G2] = assemble_constant(coords, connect, mesh, node_to_corner, refelem);

G1red  = G1(bc_data.unkn_p, bc_data.unkn_u);
G2red  = G2(bc_data.unkn_p, bc_data.unkn_v);

X(dof.u(1) -1 + bc_data.enforced_u) = bc_data.enforced_u_value(:,1);
X(dof.v(1) -1 + bc_data.enforced_v) = bc_data.enforced_v_value(:,1);
X(dof.p(1) -1 + bc_data.enforced_p) = bc_data.enforced_p_value(:,1);
X(dof.d(1) -1 + bc_data.enforced_d) = bc_data.enforced_d_value(:,1);

Error = zeros(1, maxIter);

for iter=1:maxIter
    visc = get_viscosity(X(dof.p), visc0,0);
    source = get_source(X(dof.u), X(dof.v),0);
    [K1, ~, ~, C1, ~, ~, M12_tau, K_tau, C1_tau] = assemble_iterated(connect, coords, node_to_corner, X, mesh, dof, Gamma, refelem, visc, mu, theta, dt);

    B = -mu*K0 + C1;
    d = (M12 - C1 * M12_tau)*source;
    
    % Matrices eliminations
    K1u = K1(bc_data.unkn_u, bc_data.unkn_u);
    K1v = K1(bc_data.unkn_v, bc_data.unkn_v);    
    Bred = B(bc_data.unkn_d, bc_data.unkn_d);
    
    % Load vectors
    bu = - K1(bc_data.unkn_u, bc_data.enforced_u)*bc_data.enforced_u_value(:,1) - G1(bc_data.enforced_p, bc_data.unkn_u)'*bc_data.enforced_p_value(:,1);
    bv = - K1(bc_data.unkn_v, bc_data.enforced_v)*bc_data.enforced_v_value(:,1) - G2(bc_data.enforced_p, bc_data.unkn_v)'*bc_data.enforced_p_value(:,1);
    bp = - G1(bc_data.unkn_p, bc_data.enforced_u)*bc_data.enforced_u_value(:,1) - G2(bc_data.unkn_p, bc_data.enforced_v) *bc_data.enforced_v_value(:,1);
    dred = d(bc_data.unkn_d) - B(bc_data.unkn_d, bc_data.enforced_d) * bc_data.enforced_d_value(:,1);
    
    [K, F] = build_monolythic_steady(K1u, K1v, G1red, G2red, Bred, dred, bu, bv, bp);
    
    Xred = reduce_vec(X, bc_data, dof);
    res = F - K*Xred;
    
    dXred = K\res;
    dX = build_dX(dXred, bc_data, dof, mesh);
    
    error = norm(dX([dof.u, dof.v, dof.d])); % Ignoring pressure
    X = X + dX;
    
    fprintf('Iteration %3d. Error = %.3e\n', iter, error);
    Error(iter) = error;
    if  error < tol
        break;
    end
end


loglog(Error,'+-');

post_processing_single(coords, connect, mesh, dof, corner_to_node, X)

function dX = build_dX(dXred, bc_data, dof, mesh)   
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

function Xred = reduce_vec(X, bc_data, dof)
    u_dof = bc_data.unkn_u;
    v_dof = dof.v(1) - 1 + bc_data.unkn_v;
    p_dof = dof.p(1) - 1 + bc_data.unkn_p;
    d_dof = dof.d(1) - 1 + bc_data.unkn_d;
    
    Xred = X([u_dof, v_dof, p_dof, d_dof]);
end

function [K, F] = build_monolythic_steady(K1u, K1v, G1, G2, B, d, bu, bv, bp)
    zuv = sparse(size(K1u,1), size(K1v,1));
    zpp = sparse(size(G1,1), size(G1,1));
    zud = sparse(size(K1u,1), size(B,1));
    zvd = sparse(size(K1v,1), size(B,1));
    zpd = sparse(size(G1,1), size(B,1));
    
    K = [ K1u    zuv    G1'     zud
          zuv'   K1v    G2'     zvd
          G1     G2     zpp     zpd
          zud'   zvd'   zpd'      B  ];
      
      
    F = [bu;
         bv;
         bp;
          d];
end