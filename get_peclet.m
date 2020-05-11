function [Pe_global, Pe] =  get_peclet(coords, connect, X, refelem, mu)
    
    n_elems = size(connect,2);
    n_nodes = size(coords,2);
    
    modU = sqrt(X(1:n_nodes).^2 + X(n_nodes+1:2*n_nodes).^2);
    
    Pe = zeros(n_nodes,1);
    
    for el = 1:n_elems
        nodes = connect(:,el);
        x = coords(:,nodes);
        
        jacobian = refelem.jacobian(x);
        detJ = det(jacobian.calc(jacobian,0,0));
        area = detJ * 4;
        h = sqrt(area);
        
        mean_u = mean(modU(nodes));
        
        Pe(el) = mean_u * h / (2 * mu);
        
    end
    Pe_global = max(Pe);
end