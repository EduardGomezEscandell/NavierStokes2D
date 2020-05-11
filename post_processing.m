function post_processing(coords, X_history, duration)
    n_nodes = size(coords,2);
    n_steps = size(X_history,2);
    width = max(coords(1,:));
    height = max(coords(2,:));
    
    T = delaunay(coords(1,:), coords(2,:));
    
    for step=1:n_steps
        t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), ...
               X_history(2*n_nodes+1:3*n_nodes,step));
        t.EdgeColor = 'None';
        shading interp
        hold on
        quiver(coords(1,:)',                ...
               coords(2,:)',                ...
               X_history(1:n_nodes,step),   ...
               X_history(n_nodes+1:2*n_nodes,step),'r');
        hold off
        view(2)
        axis([-1 width+1 -1 height+1]);
        title(['t=',sprintf('%.3f s',step/n_steps * duration)]);
        drawnow;
        pause(1/30);
    end
end