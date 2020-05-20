function post_processing_single(coords, connect, mesh, dof, corner_to_node, X)  
    
    width = max(coords(1,:));
    height = max(coords(2,:));
    modU = sqrt(X(dof.u).^2 + X(dof.v).^2);
    
    subplot(1,2,1);
    
    plotX = reshape(coords(1,:),[mesh.cols*2+1, mesh.rows*2+1])';
    plotY = reshape(coords(2,:),[mesh.cols*2+1, mesh.rows*2+1])';
    plotU = reshape(modU,[mesh.cols*2+1, mesh.rows*2+1])';
    t = surf(plotX, plotY, zeros(size(plotX)), plotU);  
    t.EdgeColor = 'None';
    shading interp

    hold on
    
    quiver(coords(1,:)', coords(2,:)', X(dof.u), X(dof.v),'r');
%     
%     streamX = zeros(mesh.cols, mesh.rows);
%     streamY = zeros(mesh.cols, mesh.rows);
%     streamU = zeros(mesh.cols, mesh.rows);
%     streamV = zeros(mesh.cols, mesh.rows);
%     n_lines = 10;
%     xstart = width * ones(n_lines,1);
%     ystart = linspace(0, height, n_lines);
%     s = 1;
%     for i = 1:(2*mesh.rows+1)
%         for j = 1:(2*mesh.cols+1)
%            streamX(i,j) = coords(1,s);
%            streamY(i,j) = coords(2,s);
%            streamU(i,j) = X(dof.u(1) + s - 1);
%            streamV(i,j) = X(dof.v(1) + s - 1);
%            s = s + 1;
%         end
%     end
%     
%     streamline(streamX, streamY, -streamU, -streamV, xstart, ystart);
       
    hold off
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    c=colorbar('southoutside');
    ylabel(c,'Velocity');
    title('Velocity field');
    axis equal
    
    subplot(1,2,2);
    plotX = reshape(coords(1,corner_to_node),[mesh.cols+1, mesh.rows+1])';
    plotY = reshape(coords(2,corner_to_node),[mesh.cols+1, mesh.rows+1])';
    plotP = reshape(X(dof.p),[mesh.cols+1, mesh.rows+1])';
    t = surf(plotX, plotY, zeros(size(plotX)), plotP);  
    t.EdgeColor = 'None';
    shading interp
    c=colorbar('southoutside');
    ylabel(c,'Pressure');
    view(2)
    axis([-.2 width+.2 -.2 height+.2]);
    title('Pressure');
    axis equal
    
%     subplot(1,3,3);
%     t = trisurf(T1, coords(1,corner_to_node)', coords(2,corner_to_node)', zeros(mesh.corners,1), X(dof.d));
%     t.EdgeColor = 'None';
%     shading interp
%     c=colorbar('southoutside');
%     ylabel(c,'Concentration');
%     view(2)
%     axis([-.2 width+.2 -.2 height+.2]);
%     title('Concentration');
end