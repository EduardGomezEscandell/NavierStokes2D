function post_processing_single(coords, X)
    n_nodes = size(coords,2);
    
    width = max(coords(1,:));
    height = max(coords(2,:));

    T = delaunay(coords(1,:), coords(2,:));

    t = trisurf(T, coords(1,:)', coords(2,:)', zeros(n_nodes,1), ...
           X(2*n_nodes+1:3*n_nodes));
    t.EdgeColor = 'None';
    shading interp
    hold on
    quiver(coords(1,:)',                ...
           coords(2,:)',                ...
           X(1:n_nodes),   ...
           X(n_nodes+1:2*n_nodes),'r');
    hold off
    view(2)
    axis([-1 width+1 -1 height+1]);
    colorbar;
end