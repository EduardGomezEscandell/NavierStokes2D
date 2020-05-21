function post_processing(coords, X_history, duration, dof, mesh,  corner_to_node)
    
    dt = duration / mesh.steps;

    width = max(coords(1,:));
    height = max(coords(2,:));
    
    T1 = delaunay(coords(1,corner_to_node), coords(2,corner_to_node));
    T2 = delaunay(coords(1,:), coords(2,:));
    
    time = 0;
    dt = duration / mesh.steps;
    
    modU = sqrt(X_history(dof.u,:).^2 + X_history(dof.v,:).^2);
    
    maxU = max(max(modU));
    maxp = max(max(X_history(dof.p,:)));
    maxd = max(max(X_history(dof.d,:)));
    
    minU = min(min(modU));
    minp = min(min(X_history(dof.p,:)));
    mind = min(min(X_history(dof.d,:)));
    
    for step=1:mesh.steps+1
        tic;
        X = X_history(:,step);
        modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
        
        subplot(1,3,1);

        t = trisurf(T2, coords(1,:)', coords(2,:)', zeros(mesh.nodes,1), modU);

        t.EdgeColor = 'None';
        shading interp

        hold on

        quiver(coords(1,:)', coords(2,:)', X(dof.u), X(dof.v),'r');

        hold off
        view(2)
        caxis([minU, maxU]);
        axis([-.2 width+.2 -.2 height+.2]);
        c=colorbar('southoutside');
        ylabel(c,'Velocity');
        title('Velocity field');

        subplot(1,3,2);
        t = trisurf(T1, coords(1,corner_to_node)', coords(2,corner_to_node)', zeros(mesh.corners,1), X(dof.p));
        t.EdgeColor = 'None';
        shading interp
        caxis([minp, maxp]);
        c=colorbar('southoutside');
        ylabel(c,'Pressure');
        view(2)
        axis([-.2 width+.2 -.2 height+.2]);
        title(sprintf('Pressure @ t =%8.5fs',time));

        subplot(1,3,3);
        t = trisurf(T2, coords(1,:)', coords(2,:)', zeros(mesh.nodes,1), X(dof.d));
        t.EdgeColor = 'None';
        shading interp
        c=colorbar('southoutside');
        ylabel(c,'Concentration');
        view(2)
        axis([-.2 width+.2 -.2 height+.2]);
        title('Concentration');
        
        drawnow;
        pause(max(0,dt-toc));
        time = time + dt;
    end
end