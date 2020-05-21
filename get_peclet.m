function Pe = get_peclet(coords, connect, mesh, dof, X, refelem, mu)
    
    modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
    
    Pe = zeros(mesh.nodes,1);
    elements_per_node = zeros(mesh.nodes,1);
    
    for el = 1:mesh.elems
        nodes = connect(:,el);
        x = coords(:,nodes);
        
        jacobian = refelem.Q2.jacobian(x);
        detJ = det(jacobian.calc(jacobian,0,0));
        area = detJ * 4;
        h = sqrt(area);
        
        mean_u = mean(modU(nodes));
        
        Pe_element = mean_u * h / (2 * mu);
        Pe(nodes) = Pe(nodes) + Pe_element;
        elements_per_node(nodes) = elements_per_node(nodes) + 1;
    end
    Pe = Pe ./ elements_per_node;
end