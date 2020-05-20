function data_histogram(X,variable)
    n_nodes = length(X) / 4;
    switch variable
        case 'u'
            dof = 1:n_nodes;
        case 'v'
            dof = n_nodes+1:2*n_nodes;
        case 'p'
            dof = 2*n_nodes+1:3*n_nodes;
        case 'd'
            dof = 3*n_nodes+1:4*n_nodes;
        otherwise
            error('Unrecognized variable');
    end
    bar(sort(X(dof)));
    title(['Histogram of ', variable]);
end