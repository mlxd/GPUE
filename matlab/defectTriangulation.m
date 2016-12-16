function [h,DefCount] = defectTriangulation(x, y, dx, dimSize, plotit, radius)
%Plots delaunay triangulation of input data points, with (5,7) (3,9) and
%(4,8) dislocations indicated. Can be used to calculate the defects numbers with
%no plotting by setting plotit to 0
%   radius defines the distance from 0,0, for which the defects will be
%considered. This can be used to avoid defects that arise naturally from
%the edge of the triangulation
%   h is a handle for generated plot
%   DefCount is a vector of the total number of defects counted for each
%   type with the index representing the number of connected edges

%some useful colours
RGB=[0    0.4470    0.7410
0.8500    0.3250    0.0980
0.9290    0.6940    0.1250
0.4940    0.1840    0.5560
0.4660    0.6740    0.1880
0.3010    0.7450    0.9330
0.6350    0.0780    0.1840];

%Rescale values to be centred on zero, and in meters
x=dx*(x-dimSize/2)';
y=dx*(y-dimSize/2)';

%Calculate delaunay triangulation, and setup figure
DT = delaunayTriangulation([x y]);
if plotit
    triplot(DT, 'color',[0.6 0.6 0.6],'LineWidth',2); 
    %triplot(DT, 'color',RGB(1,:),'LineWidth',3); 

    hold on;
    %set(gca,'Color',[0 0 0]);
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'DefaultTextInterpreter','Latex')
    set(gca,'FontName','Latin Modern Roman','FontSize',24);
    axis square;
    h=plot(x,y,'ko','MarkerSize',10);
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex')

    %Defect marker size
    MarkerSize = 18;
else
    h=-1;
end
%Loop over vortices and check the number of edges per index, and place
%appropriate marker.
for jj=1:9
   DefCount(1,jj)=0;
end
% idx 1 is for time, idx 2 is for type. 
for ii=1:length(x)
    for jj=1:9
        DefCount(ii+1,jj)=0;
    end
   
    if sqrt(sum([x(ii),y(ii)].^2)) < radius %%ignore edges
        if (length(DT.vertexAttachments{ii})==6)
            DefCount(ii+1,6)= DefCount(ii+1,6) +1;
        end

        if (length(DT.vertexAttachments{ii})==5)
            if plotit
                plot(x(ii),y(ii),'p','MarkerEdgeColor','k','MarkerSize',MarkerSize-2,'MarkerFaceColor','w','LineWidth',1.5);
            end
            DefCount(ii+1,5)= DefCount(ii+1,5) +1;
        elseif length(DT.vertexAttachments{ii})==7
            if plotit
                plot(x(ii),y(ii),'h','MarkerEdgeColor','k','MarkerSize',MarkerSize,'MarkerFaceColor',[0.3 0.3 0.3],'LineWidth',1.5);
            end
            DefCount(ii+1,7)= DefCount(ii+1,7) +1;
        end

        if (length(DT.vertexAttachments{ii})==3)
            if plotit
                plot(x(ii),y(ii),'^','MarkerEdgeColor','k','MarkerSize',MarkerSize-2,'MarkerFaceColor',RGB(4,:),'LineWidth',1.5);
            end
            DefCount(ii+1,3)= DefCount(ii+1,3) +1;
        elseif length(DT.vertexAttachments{ii})==9
            if plotit
                plot(x(ii),y(ii),'v','MarkerEdgeColor','k','MarkerSize',MarkerSize-2,'MarkerFaceColor',RGB(5,:), 'LineWidth',1.5);
            end
            DefCount(ii+1,9)= DefCount(ii+1,9) +1;
        end

        if (length(DT.vertexAttachments{ii})==4)
            if plotit
                plot(x(ii),y(ii),'d','MarkerEdgeColor','k','MarkerSize',MarkerSize-2,'MarkerFaceColor',RGB(2,:),'LineWidth',1.5);
            end
            DefCount(ii+1,4)= DefCount(ii+1,4) +1;

        elseif length(DT.vertexAttachments{ii})==8
            if plotit
                plot(x(ii),y(ii),'o','MarkerEdgeColor','k','MarkerSize',MarkerSize-2,'MarkerFaceColor',RGB(6,:),'LineWidth',1.5);
            end
            DefCount(ii+1,8)= DefCount(ii+1,8) +1;
        end
    end
end
DefCount = sum(DefCount);

end