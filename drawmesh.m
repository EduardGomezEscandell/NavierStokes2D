function drawmesh(coords, connect)
    n_elems = length(connect);
    n_nodes = length(coords);
    
    for el = 1:n_elems
        X = coords(1,connect(1:4,el));
        Y = coords(2,connect(1:4,el));
        X(5) = X(1);
        Y(5) = Y(1);
        plot(X,Y,'k');
        hold on
        text(mean(X),mean(Y),num2str(el),'Color','k');
    end
    
    scatter(coords(1,:),coords(2,:),'square','filled','b')
    
    labels = cell(n_nodes,1);
    for i = 1:n_nodes
        labels{i} = sprintf(' %d',i);
    end
    
    text(coords(1,:),coords(2,:),labels,'Color','b');
    hold off
end