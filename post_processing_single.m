function post_processing_single(coords, connect, X)
    n_nodes = size(coords,2);
    
    u_dof = 1:n_nodes;
    v_dof = n_nodes+1:2*n_nodes;
    p_dof = 2*n_nodes+1:3*n_nodes;
    d_dof = 3*n_nodes+1:4*n_nodes;
    
    width = max(coords(1,:));
    height = max(coords(2,:));
    modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
    
    T = delaunay(coords(1,:), coords(2,:));
    
    subplot(1,3,1);

    t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), modU);
    
    t.EdgeColor = 'None';
    shading interp

    hold on
    
    quiver(coords(1,:)', coords(2,:)', X(u_dof), X(v_dof),'r');
       
    hold off
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    c=colorbar('southoutside');
    ylabel(c,'Velocity');
    title('Velocity field');
    
    subplot(1,3,2);
    t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), X(p_dof));
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Pressure');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Pressure');
    
    subplot(1,3,3);
    t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), X(d_dof));
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Concentration');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Concentration');
end