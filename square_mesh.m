function [coords,connect, corner_to_node, node_to_corner] = square_mesh(width, height, x_elems, y_elems, degree)
    % Returns the coordinates and connectivities of a square mesh on domain
    % [0,width]x[0,height]
    
    nodes_per_row = x_elems*degree + 1;
    nodes_per_col = y_elems*degree + 1;
    n_nodes = nodes_per_col * nodes_per_row;
    
    X = linspace(0, width, nodes_per_row);
    Y = linspace(0, height, nodes_per_col);
    
    coords = zeros(2, n_nodes);
    
    corner_to_node = [];
    node_to_corner = zeros(1, n_nodes);
    
    k = 1;
    for row=1:nodes_per_col
        for col=1:nodes_per_row
            node_to_corner(k) = -1;
            if mod(row,2) == 1 && mod(col,2)==1
                corner_to_node(end+1) = k;
                node_to_corner(k) = length(corner_to_node);
            end
            
            coords(:,k) = [X(col), Y(row)]';
            k = k + 1;
        end
    end
    
    switch degree
        case 1
             % Local coordinates:
            % 4 - 3
            % |   |
            % 1 - 2
            bottom_lefters = 1:n_nodes;
            
            must_stay = logical( (mod(bottom_lefters,nodes_per_row) ~= 0) ...
                               .*(bottom_lefters < (n_nodes - nodes_per_row)));
            
            bottom_lefters = bottom_lefters(must_stay);


            connect = [bottom_lefters;                % Bottom left
                       bottom_lefters+1;              % Bottom right
                       bottom_lefters+nodes_per_row+1 % Top right
                       bottom_lefters+nodes_per_row]; % Top left
        case 2
            % Local coordinates:
            % 4 - 7 - 3
            % |       |
            % 8   9   6
            % |       |
            % 1 - 5 - 2
            bottom_lefters = 1:n_nodes;
            row = ceil(bottom_lefters / nodes_per_row);
            col = mod(bottom_lefters-1, nodes_per_row)+1;
            
            must_stay = logical( (mod(row,2)==1) ...
                               .*(mod(col,2)==1) ...
                               .*(col ~= nodes_per_row) ...
                               .*(row ~= nodes_per_col));
             
            bottom_lefters = bottom_lefters(must_stay);

            connect = [bottom_lefters;                    % Bottom left
                       bottom_lefters+2;                  % Bottom right
                       bottom_lefters+2*nodes_per_row+2   % Top right
                       bottom_lefters+2*nodes_per_row;    % Top left
                       bottom_lefters + 1;                % Bottom
                       bottom_lefters+  nodes_per_row+2;  % Right
                       bottom_lefters+2*nodes_per_row+1;  % Top
                       bottom_lefters+  nodes_per_row;    % Left
                       bottom_lefters+  nodes_per_row+1]; % Center
            
        otherwise
            error('Only first and second order implemented');
    end
end