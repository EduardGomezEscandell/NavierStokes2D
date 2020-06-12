function post_processing_single(coords, connect, mesh, dof, corner_to_node, X)  
    
    width = max(coords(1,:));
    height = max(coords(2,:));    
    
    %% Velocity
    
    plotX = reshape(coords(1,:),[mesh.cols*2+1, mesh.rows*2+1])';
    plotY = reshape(coords(2,:),[mesh.cols*2+1, mesh.rows*2+1])';
    
    modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
    plotU = reshape(modU,[mesh.cols*2+1, mesh.rows*2+1])';
    
    subplot(1,3,1);
    surf(plotX, plotY, zeros(size(plotX)), plotU);  
    shading interp

    hold on    
    quiver(coords(1,:)', coords(2,:)', X(dof.u), X(dof.v),'r');
    hold off
    view(2)
    c=colorbar('southoutside');
    ylabel(c,'Velocity');
    title('Velocity field');
    axis tight
    
    %% Pressure
    subplot(1,3,2);
    plotX = reshape(coords(1,corner_to_node),[mesh.cols+1, mesh.rows+1])';
    plotY = reshape(coords(2,corner_to_node),[mesh.cols+1, mesh.rows+1])';
    plotP = reshape(X(dof.p),[mesh.cols+1, mesh.rows+1])';
    surf(plotX, plotY, plotP, plotP);
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Pressure');
    view(2)
    caxis([prctile(X(dof.p),10), prctile(X(dof.p),90)]);
    title('Pressure');
    
    %% Concentration
    subplot(1,3,3);
    plotD = reshape(X(dof.d),[mesh.cols+1, mesh.rows+1])';
    t = surf(plotX, plotY, plotD, plotD);
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Concentration');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Concentration');
    axis tight;
end