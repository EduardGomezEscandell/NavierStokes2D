function refelem = set_reference_element(degree)
    switch degree
        case 1
            refelem.order = 1;
            refelem.shape = 'Q';
            refelem.jacobian = @jacobian4;
            
            N{1} = @(x,y)((1-x)*(1-y)/4);
            N{2} = @(x,y)((1+x)*(1-y)/4);
            N{3} = @(x,y)((1+x)*(1+y)/4);
            N{4} = @(x,y)((1-x)*(1+y)/4);
            
            gradN{1} = @(x,y)([ (y-1)/4,  (x-1)/4]');
            gradN{2} = @(x,y)([ (1-y)/4, -(1+x)/4]');
            gradN{3} = @(x,y)([ (1+y)/4,  (1+x)/4]');
            gradN{4} = @(x,y)([-(1+y)/4,  (1-x)/4]');
            
            gauss_p = [1,-1] / sqrt(3);
            refelem.gauss_p = gauss_p;
            refelem.weights = ones(2);
            
            p = 1;
            for xi = gauss_p                
                for eta=gauss_p
                    for i=1:length(N)
                        refelem.N(i,p) = N{i}(xi,eta);
                        refelem.gradN(:,i,p) = gradN{i}(xi,eta);
                    end
                    p = p+1;
                end
            end
            
        case 2
            refelem.order = 2;
            refelem.shape = 'Q';
            refelem.jacobian = @jacobian9;
            
            % Shape functions
            N{1} = @(x,y) ( 1/4 * (1-x) * (1-y) * x*y);
            N{2} = @(x,y) (-1/4 * (1+x) * (1-y) * x*y);
            N{3} = @(x,y) ( 1/4 * (1+x) * (1+y) * x*y);
            N{4} = @(x,y) (-1/4 * (1-x) * (1+y) * x*y);
            N{5} = @(x,y) (-1/2 * (1-x^2) * (1-y) * y);
            N{6} = @(x,y) ( 1/2 * (1+x) * (1-y^2) * x);
            N{7} = @(x,y) ( 1/2 * (1-x^2) * (1+y) * y);
            N{8} = @(x,y) (-1/2 * (1-x) * (1-y^2) * x);
            N{9} = @(x,y) ((1-x^2) * (1-y^2));
            
            % Gradients of shape functions
            gradN{1} = @(x,y) (1/4*[y*(1-2*x)+y^2*(2*x-1), x*(1-2*y)+x^2*(2*y-1)]');
            gradN{2} = @(x,y) (-1/4*[y*(1+2*x)-y^2*(2*x+1), x*(1-2*y)+x^2*(1-2*y)]');
            gradN{3} = @(x,y) (1/4*[(1+2*x)*(y+y^2), (1+2*y)*(x+x^2)]');
            gradN{4} = @(x,y) (-1/4*[(1-2*x)*(y+y^2), (1+2*y)*(x-x^2)]');
            gradN{5} = @(x,y) (-1/2*[2*x*y*(y-1), 1-2*y-x^2+2*x^2*y]');
            gradN{6} = @(x,y) (1/2*[1-y^2+2*x-2*x*y^2, -2*x*y*(x+1)]');
            gradN{7} = @(x,y) (1/2*[-2*x*y*(y+1), 1+2*y-x^2-2*x^2*y]');
            gradN{8} = @(x,y) (-1/2*[1-y^2-2*x+2*x*y^2, 2*x*y*(x-1)]');
            gradN{9} = @(x,y) ([2*x*(y^2-1), 2*y*(x^2-1)]');
            
            gauss_p = [1,-1] / sqrt(3);
            refelem.gauss_p = gauss_p;
            refelem.weights = ones(2);
            
            p = 1;
            for xi = gauss_p                
                for eta=gauss_p
                    
                    for i=1:length(N)
                        refelem.N(i,p) = N{i}(xi,eta);
                        refelem.gradN(:,i,p) = gradN{i}(xi,eta);
                    end
                    p = p+1;
                end
            end
            
            
        otherwise
            error('Higher-than-quadratic shape functions unavailable');
    end
end

function jacobian = jacobian4(X)
    j0(1,:) = (- X(:,1) + X(:,2) + X(:,3) - X(:,4))';
    j0(2,:) = (- X(:,1) - X(:,2) + X(:,3) + X(:,4))';
    j1(1,:) = (+ X(:,1) - X(:,2) + X(:,3) - X(:,4))';
    
    jacobian.j0 = j0;
    jacobian.j1 = j1;
    jacobian.calc = @(jacobian, xi,eta)(0.25 * jacobian.j0 + [xi;eta]*jacobian.j1);
end

function jacobian = jacobian9(X)
    
end