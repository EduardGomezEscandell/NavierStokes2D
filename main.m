clear all;clc;

% Domain
width = 2;
height = 1;
duration = 5;

% Discretization
x_elems = 4;
y_elems = 4;
n_steps = 50;

degree = 1; %(space)

% Non-linearity
maxIter = 50;
tol = 1e-5;

% Physical properties
visc0 = 1;
mu = 1;

% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
% drawmesh(coords, connect)
refelem = set_reference_element(degree);

n_nodes = length(coords);
n_elems = length(connect);
nodes_per_elem = size(connect,1);

% Solution vectors
X_history = zeros(n_nodes*3, n_steps);
X = X_history(:,1);

% First iteeration
source = zeros(n_nodes,1);

% Initializing
dt = duration / n_steps;
K_global = sparse(n_nodes*3, n_nodes*3); % ux, uy, p, r
F_global = sparse(n_nodes*3,1);

% Assembling

for step = 2:n_steps
        
    for iter=1:maxIter
        
        modU = sqrt(X(1:n_nodes).^2 + X(n_nodes+1:n_nodes*2).^2);
        source_prev = source;
        source = 1./(1 + exp(-10*(modU - 0.5)));
        rho = X(2*n_nodes+1:3*n_nodes);
        visc = visc0 + visc0 ./ (1 + exp(-10*(rho - 0.5)));
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            a =  [X(nodes), X(nodes+n_nodes)];
            s = source(nodes,:);
            s_prev = source_prev(nodes,:);
            v = visc(nodes,:);
            [K, K_v, C, M, F, F_prev] = FEM_matrices(local_coordinates,refelem, a, s, s_prev, v);
            
            % Assembly of Solvent
            % A*u = 0
            A =  [dt*M + 0.5*K_v, zeros(size(K)); zeros(size(K)), dt*M + 0.5*K_v];
            b = -[K_v, zeros(size(K)); zeros(size(K)), K_v] * X_history([nodes; nodes+n_nodes], step-1);
            
            % Assembly of solute
            % B*u = S
            B = (dt*M + 0.5*C + mu*K);
            S = -(C + K)*rho(nodes,:) + M*(F_prev + F)/2; % Assembly of Source term
            
            % Joining two systems
            % K*u = F
            K_local = [                                     A, zeros(2*nodes_per_elem,nodes_per_elem);
                       zeros(nodes_per_elem,2*nodes_per_elem),                                      B];
                   
            F_local = [b; S];
                      
            % Assembling
            nodes = [nodes;
                     nodes + n_nodes;
                     nodes + n_nodes*2];

            K_global(nodes, nodes) = K_global(nodes, nodes) + K_local;
            F_global(nodes) = F_global(nodes) + F_local; 
        end
        
        % Newton-Rhapson
        R = K_global*X - F_global;
        
        % Calculating jacobian
        
        % Newton Rhapson
        X_new = X - dRdX \ R;
        
        % Convergence condition
        error = norm(R);
        if error < tol
            break;
        end
        X = X_new;
    end
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    Error(step) = error;
    
    if(mod(step - 1,10) == 0)
        fprintf('Step %3d of %3d completed\n',step, n_steps);
    end
    
    S_prev = S;
    X_history(:,step) = X_new;
end

post_porcessing(coords, X_history, duration)