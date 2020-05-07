% Domain
width = 6;
height = 2;
duration = 5;

% Discretization
x_elems = 12;
y_elems = 4;
n_steps = 10;

degree = 1; %(space)

% Non-linearity
maxIter = 50;
tol = 1e-8;

% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
drawmesh(coords, connect)
refelem = set_reference_element(degree);

n_nodes = length(coords);
n_elems = length(connect);
nodes_per_elem = size(connect,1);

% Solution vectors
X_history = zeros(n_nodes*4, n_steps);
X = X_history(:,1);

K_global = sparse(n_nodes*4, n_nodes*4); % ux, uy, p, r
F_global = sparse(n_nodes*4,1);

% Assembling
tic;
for step = 2:n_steps
        
    for iter=1:maxIter
            
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            k = [1,0;0,1];
            c = [1;0];
            m = 1;
            
            [K, C, M] = FEM_matrices(local_coordinates,refelem, k,c,m);
            
            % Assembly
            A = ones(nodes_per_elem*3,nodes_per_elem*3); % Assembly of A
            B = ones(nodes_per_elem,nodes_per_elem);   % Assembly of B
            S = ones(nodes_per_elem,1);   % Assembly of Source term
            
            K_local = [                                     A, zeros(3*nodes_per_elem,nodes_per_elem);
                       zeros(nodes_per_elem,3*nodes_per_elem),                                      B];
                   
            F_local = [zeros(12, 1); S];
                      
            % Assembling
            nodes = [nodes;
                     nodes + n_nodes;
                     nodes + n_nodes*2;
                     nodes + n_nodes*3];

            K_global(nodes, nodes) = K_global(nodes, nodes) + K_local;
            F_global(nodes) = F_global(nodes) + F_local; 
        end
        
        % Newton-Rhapson
        G = K_global*X - F_global;
        J = ones(size(K_global)); % Calculation of Jacobian
        X_new = X - J*G;
        
        % Convergence condition
        if norm(X_new - X) < tol
            break;
        end
        X = X_new;
    end
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    X_history(:,step) = X_new;
end
toc
% post_porcessing(coords, X_history, duration)