function drawmesh(coords, connect, draw_labels, varargin)
    
    if nargin == 2
        draw_labels = true;
    elseif nargin < 2 || nargin > 3
        error('Wrong number of arguments');
    end
    
    n_elems = length(connect);
    n_nodes = length(coords);
    
    for el = 1:n_elems
        X = coords(1,connect(1:4,el));
        Y = coords(2,connect(1:4,el));
        X(5) = X(1);
        Y(5) = Y(1);
        plot(X,Y,'k');
        hold on
        if draw_labels
            text(mean(X),mean(Y),num2str(el),'Color','k');
        end
    end
    
    if draw_labels
        scatter(coords(1,:),coords(2,:),'square','filled','b')
        
        labels = cell(n_nodes,1);
        for i = 1:n_nodes
            labels{i} = sprintf(' %d',i);
        end
    
        text(coords(1,:),coords(2,:),labels,'Color','b');
    end
    hold off
end