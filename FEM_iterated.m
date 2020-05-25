function local_mat = FEM_iterated(local_coords, vel, conc, visc, mu, theta, dt, linear_elem, quadra_elem)
    % First Order
    C1 = zeros(4,4);
    C21 = zeros(4,4);
    C22 = zeros(4,4);
    
    % Second Order
    K1 = zeros(9,9);
    K21 = zeros(9,9);
    K22 = zeros(9,9);
    Q = zeros(9,9);
    
    % Hybrid
    K_tau = zeros(4,4);
    C1_tau = zeros(4,4);
    M12_tau = zeros(4,9);
    
    jacobian1 = linear_elem.jacobian(local_coords);
    jacobian2 = quadra_elem.jacobian(local_coords);
    
    detJ = det(jacobian1.calc(jacobian1,0,0));
    area = 4*detJ;
    h = sqrt(area);
       
    % Integration
    p = 1;
    for gp1 = linear_elem.gauss
        xi = gp1(1);
        for gp2 = linear_elem.gauss
            eta = gp2(1);
            w = gp1(2) * gp2(2);

            % Linear Shape functions
            J1 = jacobian1.calc(jacobian2,xi,eta);
            gradN1 = J1\linear_elem.gradN(:,:,p);
            N1 = linear_elem.N(:,p)';
            
            % Quadratic Shape functions
            J2 = jacobian2.calc(jacobian2,xi,eta);
            gradN2 = J2\quadra_elem.gradN(:,:,p);
            N2 = quadra_elem.N(:,p)';

            % Gradients
            grad_u = gradN2 * vel(:,1);
            grad_v = gradN2 * vel(:,2); 
            grad_dRho= gradN1 * conc;

            % Diffusion
            K1  =  K1 + w * gradN2' * (N1*visc) * gradN2; % K1(nu)
            K21 = K21 + w * gradN2' * grad_u * N2;           % K2(u)
            K22 = K22 + w * gradN2' * grad_v * N2;           % K2(v)

            % Convection
            C1  =  C1 + w * N1' * ((N2*vel * gradN1)); 
            C21 = C21 + w * (N1' * N1) * grad_dRho(1);      % C21(rho)
            C22 = C22 + w * (N1' * N1) * grad_dRho(2);      % C22(rho)
            
            % Stabilization in concentration
            a = N2*vel;
            tau = (1/(theta*dt) + 2*norm(a)/h + 4*mu/h^2)^-1;
%             tau = (1/(theta*dt)^2 + (2*norm(a)/h)^2 + 9*(4*mu/h^2)^2)^(-1/2);
            
            M12_tau= M12_tau + w * tau * (N1'*N2);
            K_tau  = K_tau   + w * tau * (gradN1' * gradN1);
            C1_tau = C1_tau  + w * tau * N1' * ((N2*vel * gradN1)); 

            p = p+1;
        end
    end

    local_mat.K1  = K1  * detJ;
    local_mat.K21 = K21 * detJ;
    local_mat.K22 = K22 * detJ;
    
    local_mat.C1  = C1  * detJ;
    local_mat.C21 = C21 * detJ;
    local_mat.C22 = C22 * detJ;
    
    local_mat.M12_tau = M12_tau * detJ;
    local_mat.K_tau   = K_tau   * detJ;
    local_mat.C1_tau  = C1_tau  * detJ;
        
end