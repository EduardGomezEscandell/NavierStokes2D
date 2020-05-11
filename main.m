clear variables;
clc;

% Domain
width = 2;
height = 3;
duration = 1;

% Discretization
x_elems = 8;
y_elems = 12;
degree = 1; %(space)

n_steps = 10;
theta = 1/2;

method = 1;
penalty = 1e5;
% 1: Picard
% 2: Newton-Raphson

% Non-linearity
maxIter = 50;
tol = 1e-5;

% Physical properties
visc0 = 100;
mu = 100;
omega = 1;

% Generating mesh
[coords,connect] = square_mesh(width, height, x_elems, y_elems, degree);
refelem = set_reference_element(degree);

% Shorthand notation
n_nodes = length(coords);
n_elems = length(connect);

u_dof = 1:n_nodes;
v_dof = n_nodes+1:2*n_nodes;
p_dof = 2*n_nodes+1:3*n_nodes;
d_dof = 3*n_nodes+1:4*n_nodes;

dt = duration / n_steps;

% Solution vectors
X_history = zeros(n_nodes*4, n_steps);

% Selecting boundary conditions
[free_dof, X_history, H, e] =  boundary_conditions(coords, X_history, duration, omega);

% Assembling

for step = 2:n_steps
    X = X_history(:,step-1);
    
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    source = 1./(1 + exp(-10*(modU - 0.5)));
    visc = visc0 + visc0 ./ (1 + exp(-10*(X(d_dof) - 0.5)));
        
    % Obtaining arrays and matrices that need not iteration
    b1 = zeros(n_nodes,1);
    b2 = zeros(n_nodes,1);
    d_non_iter = zeros(n_nodes,1);
    
    for el=1:n_elems
        nodes = connect(:,el);
        local_coordinates = coords(:, nodes);
        
        x  = [X(nodes), X(nodes+n_nodes), X(nodes+2*n_nodes), X(nodes+3*n_nodes)];
        v = visc(nodes);
        s = source(nodes);
        
        [M, K, K1, C1, Q1, Q2] = non_iterated_arrays(local_coordinates, x, v, refelem);
        
        b = (M + dt*(1-theta)*K1);
        d = (M - dt*(1-theta)*(C1 + mu*K));
        
        b1(nodes) = b1(nodes) + b * x(:,1) + dt * (1-theta)*Q1*x(:,3);
        b2(nodes) = b2(nodes) + b * x(:,2) + dt * (1-theta)*Q1*x(:,3);
        
        d_non_iter(nodes)  =  d_non_iter(nodes) ...
                  + d*x(:,4) + dt * (1-theta)*M*s;
    end
    
    X = X_history(:,step);
    X(free_dof) = X_history(free_dof,step-1);
    
    % Iterating
    for iter=1:maxIter
        M   = sparse(n_nodes,n_nodes);
        
        K   = sparse(n_nodes,n_nodes);
        K1  = sparse(n_nodes,n_nodes);
        K21 = sparse(n_nodes,n_nodes);
        K22 = sparse(n_nodes,n_nodes);
        
        C1  = sparse(n_nodes,n_nodes);
        C21 = sparse(n_nodes,n_nodes);
        C22 = sparse(n_nodes,n_nodes);
        
        Q1  = sparse(n_nodes,n_nodes);
        Q2  = sparse(n_nodes,n_nodes);
        
        modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
        source = 1./(1 + exp(-10*(modU - 0.5)));
        visc = visc0 + visc0 ./ (1 + exp(-10*(X(d_dof) - 0.5)));
        
        for el=1:n_elems
            nodes = connect(:,el);
            local_coordinates = coords(:, nodes);
            
            x  = [ X(nodes), X(nodes+n_nodes), X(nodes+2*nodes)];
            s = source(nodes);
            v = visc(nodes);
                  
            local_mat = FEM_matrices(local_coordinates, refelem, x, v);
            
            % Mass matrix
            M(nodes,nodes) 	 = M(nodes,nodes)   + local_mat.M;
            
            % Difusion matrices
            K(nodes,nodes) 	 = K(nodes,nodes)   + local_mat.K;
            K1(nodes,nodes)  = K1(nodes,nodes)  + local_mat.K1;
            K21(nodes,nodes) = K21(nodes,nodes) + local_mat.K21;
            K22(nodes,nodes) = K22(nodes,nodes) + local_mat.K22;
            
            % Convection matrices
            C1(nodes,nodes)  = C1(nodes,nodes)  + local_mat.C1;
            C21(nodes,nodes) = C21(nodes,nodes) + local_mat.C21;
            C22(nodes,nodes) = C22(nodes,nodes) + local_mat.C22;
            Q1(nodes,nodes)  = Q1(nodes,nodes)  + local_mat.Q1;
            Q2(nodes,nodes)  = Q2(nodes,nodes)  + local_mat.Q2;
        end
        
        A = M + dt*theta*K1;
        T1 = -dt*theta*Q1;
        T2 = -dt*theta*Q2;
        
        B = M + dt*theta*(C1 + mu*K);
        d = d_non_iter + dt * theta * M*source;
        
        if method==1
            % Picard
            Z = zeros(size(A));
            z = zeros(size(d));

            K = [Q1  Q2   Z  Z;
                  A   Z  T1'  Z;
                  Z   A  T2'  Z;
                  Z   Z   Z  B];

            F = [ z; b1; b2; d];
            
            % Penalty method
            K = K + penalty * (H'*H);
            F = F - penalty * H'*e(:,step);
            
%             K = K(free_dof,free_dof);
%             F = F(free_dof);
            
            R = K*X - F;
            error = norm(R);
            
            %
            post_processing_single(coords, X);
            title(sprintf('Error = %g',error)); 
            %
            
            if error < tol
                break;
            end
            
            X = K\F;
            
        else
            % Newton-Rhapson

            R(u_dof,:) = Q1*X(u_dof) + Q2*X(v_dof);
            R(v_dof,:) =  A*X(u_dof) + T1*X(p_dof) - b2;
            R(p_dof,:) =  A*X(v_dof) + T2*X(p_dof) - d;
            R(d_dof,:) =  B*X(d_dof) - d;

            R = R(free_dof);

            % Convergence condition
            error = norm(R);
            
            %
            post_processing_single(coords, X);
            title(sprintf('Error = %g',error)); 
            %
            
            if error < tol
                break;
            end

            % Calculating jacobian
            J = calc_jacobian(coords, visc0, theta, dt, X, K21, K22, C21, C22, Q1, Q2, T1, T2, M, A, B);
            J = J(free_dof, free_dof);

            % Newton Raphson
            X(free_dof) = X(free_dof) - (J \ R);
        end
    end
    
    if iter==maxIter
        warning('Maximum number of iterations reached before convergence');
    end
    
    if(mod(step - 1,10) == 0)
        fprintf('Step %3d of %3d completed\n',step, n_steps);
    end
    
    X_history(:,step) = X;
end

post_processing(coords, X_history, duration)