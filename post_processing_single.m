function post_processing_single(coords, connect, mesh, dof, corner_to_node, X)  
    
    width = max(coords(1,:));
    height = max(coords(2,:));
    modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
    
    T1 = delaunay(coords(1,corner_to_node), coords(2,corner_to_node));
    T2 = delaunay(coords(1,:), coords(2,:));
    
    subplot(1,3,1);

    t = trisurf(T2, coords(1,:)', coords(2,:)', zeros(mesh.nodes,1), modU);
    
    t.EdgeColor = 'None';
    shading interp

    hold on
    
    quiver(coords(1,:)', coords(2,:)', X(dof.u), X(dof.v),'r');
       
    hold off
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    c=colorbar('southoutside');
    ylabel(c,'Velocity');
    title('Velocity field');
    
    subplot(1,3,2);
    t = trisurf(T1, coords(1,corner_to_node)', coords(2,corner_to_node)', zeros(mesh.corners,1), X(dof.p));
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Pressure');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Pressure');
    
    subplot(1,3,3);
    t = trisurf(T1, coords(1,corner_to_node)', coords(2,corner_to_node)', zeros(mesh.corners,1), X(dof.d));
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Concentration');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Concentration');
end