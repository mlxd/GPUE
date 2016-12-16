function [area,avg_area,num_edges,var]=voronoi_cells(x,y,q,dx)
[v,c]=voronoin([x y]);
area=zeros(length(c),1);
num_edges=zeros(length(c),1);
avg_area = 0;
clf;
for jj = 1:length(c)
    if (all(c{jj}~=1) && all(sqrt((x(jj)).^2 + (y(jj)).^2)<230*dx)); % If at least one of the indices is 1,
        % then it is an open region and we can't
        % patch that.
        area(jj)=polyarea(v(c{jj},1),v(c{jj},2));
        avg_area = avg_area + area(jj);
        h = patch(v(c{jj},1),v(c{jj},2),area(jj)); % use area as color
        if q(jj)< 0
            h.EdgeColor = 'red';
            h.LineWidth = '2';
        elseif q(jj) == 0
            h.EdgeColor = 'black';
            h.LineStyle = '-'
        else
            h.EdgeColor = 'white';
            h.LineStyle = '-'
        end
        num_edges(jj) = length(c{jj});
        get(h); 
    end
end
%axis square;axis([200 800 200 800]*dx);caxis([740 860]*dx*dx);;colorbar;
hold on; plot(x,y,'r*');plot(x,y,'bo');
set(gca,'FontName','Latin Modern Roman','FontSize',22);
set(gca,'TickLabelInterpreter', 'latex');axis off;
drawnow;

avg_area=avg_area./nnz(area);
var = sum((area(area>1) - avg_area).^2);
var = var./nnz(area);