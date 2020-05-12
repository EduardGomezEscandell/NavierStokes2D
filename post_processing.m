function post_processing(coords, X_history, duration)
    n_nodes = size(coords,2);
    n_steps = size(X_history,2);
    width = max(coords(1,:));
    height = max(coords(2,:));
    
    u_dof = 1:n_nodes;
    v_dof = n_nodes+1:2*n_nodes;
    p_dof = 2*n_nodes+1:3*n_nodes;
    d_dof = 3*n_nodes+1:4*n_nodes;
    
    T = delaunay(coords(1,:), coords(2,:));
    time = 0;
    dt = duration / n_steps;
    
    modU = sqrt(X_history(u_dof,:).^2 + X_history(v_dof,:).^2);
    
    maxU = max(max(modU));
    maxp = max(max(X_history(p_dof,:)));
    maxd = max(max(X_history(d_dof,:)));
    
    minU = min(min(modU));
    minp = min(min(X_history(p_dof,:)));
    mind = min(min(X_history(d_dof,:)));
    
    for step=1:n_steps
        tic;
        X = X_history(:,step);
        modU = sqrt(X(u_dof).^2 + X(v_dof).^2);
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
        caxis([minU, maxU]);
        title(sprintf('Velocity field @ t=%g',time));

        subplot(1,3,2);
        t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), X(p_dof));
        t.EdgeColor = 'None';
        shading interp
        c=colorbar('southoutside');
        ylabel(c,'Pressure');
        view(2)
        caxis([minp, maxp]);
        axis([-.2 width+.2 -.2 height+.2]);
        title('Pressure');

        subplot(1,3,3);
        t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), X(d_dof));
        t.EdgeColor = 'None';
        shading interp
        c=colorbar('southoutside');
        ylabel(c,'Concentration');
        view(2)
        caxis([mind, maxd]);
        axis([-.2 width+.2 -.2 height+.2]);
        title('Concentration');
        drawnow;
        pause(max(0,1/30-toc));
        time = time + dt;
    end
end