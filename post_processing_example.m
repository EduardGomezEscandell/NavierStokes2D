% Domain
width = 6;
height = 2;
duration = 10;

% Discretization
x_elems = 12;
y_elems = 8;
n_steps = 100;

degree = 1; %(space)

% Non-linearity
maxIter = 50;
tol = 1e-8;

% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
n_nodes = length(coords);

% Making up results
X_history = zeros(4*n_nodes, n_steps);

for step=1:n_steps
    X_history(          1:  n_nodes,step) = ones(n_nodes,1);
    X_history(  n_nodes+1:2*n_nodes,step) = sin(2*step/n_steps*2*pi - 0.3 * coords(1,:));
    X_history(2*n_nodes+1:3*n_nodes,step) = sin(2*step/n_steps*2*pi - 0.3 * coords(1,:));
end

% Post-processing
post_porcessing(coords, X_history, duration)
