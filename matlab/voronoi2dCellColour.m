function [p6cp6, area, avg_area, num_edges, var]=voronoi2dCellColour(x,y,radius,edgeCol,dx,colScheme)
%vorCellColour Determine the Voronoi diagram of the input data, and colour the
%cells with the value of the orientational correlations defined at each
%site (col=1), or with the area of the cell (col=0).
%x,y: are the locations of the sites
%radius: is the area over which to perform the triangulation. This should
%     avoid too large values to ensure the traingulation does not run off to
%     infinity
%edgeCol: Defines the Voronoi diagram edge colours for edge indexed pair. edgeCol<0 is 
%     red, edgeCol==0 is black, edgeCol>0 is white
%dx: this is the increment size of the data
%colScheme: defines the colouring scheme, 1 for orientational correlations on the site 
%     with 6-fold symmetry, 0 for cell area 
% Returns:  p6cp6 = local values of correlations
%           area = cell areas
%           avg_area = average area of cells
%           num_edges = cell edges count
%           var = variance
% 
%Testcase: Generate a 2D grid centred on 0 and calculate voronoi diagram
% x=linspace(-1,1,20); y=x;
% radius = 0.5; dx = 1; q = ones(length(x)*length(y),1);
% voronoi2dCellColour(kron(x',ones(length(y),1)),kron(ones(length(x),1),y'),0.75,zeros(length(x)^2,1),1,0);

DT = delaunayTriangulation([x y]);
[v,c]=voronoiDiagram(DT);
if colScheme==1
    p6cp6=zeros(length(c),1);
end
num_edges=zeros(length(c),1);
avg_area = 0;
area=zeros(length(c),1);

clf;
R = DT.Points;

%Credit to MathWorks support team for neighboringVertices.m
vTriAtt = vertexAttachments(DT);
for ii = 1:size(R,1)
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    verticesOfTI         = DT.ConnectivityList(vTriAtt{ii},:);
    % 3. Find all the unique vertices and remove the current vertex
    neighboursOfInternal{ii} = setdiff(unique(verticesOfTI), ii);
end

% Calculate the Voronoi diagram
for jj = 1:length(c)
    %Credit to MathWorks team for voronoi examples. https://uk.mathworks.com/help/matlab/ref/voronoi.html
    if (all(c{jj}~=1) && all(sqrt((x(jj)).^2 + (y(jj)).^2) < radius)); % If at least one of the indices is 1,
        % then it is an open region and we can't patch that.

        area(jj)=polyarea(v(c{jj},1),v(c{jj},2));
        avg_area = avg_area + area(jj);
        if colScheme==1
            [p6(jj),nn(jj)] = psi6_DT([x(jj),y(jj)],x(neighboursOfInternal{jj}),y(neighboursOfInternal{jj}));
            p6cp6(jj) = sum(conj(p6(jj)).*p6(jj));
            h = patch(v(c{jj},1),v(c{jj},2),p6cp6(jj)); % use g6 as color
        else
            h = patch(v(c{jj},1),v(c{jj},2),area(jj)); % use area as color
        end
        if edgeCol(jj)< 0
            h.EdgeColor = 'red';
            h.LineWidth = '2';
            h.LineStyle = '-';
        elseif edgeCol(jj) == 0
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

avg_area=avg_area./nnz(area);
var = sum((area(area>1) - avg_area).^2);
var = var./nnz(area);

hold on; plot(x,y,'r*');plot(x,y,'bo');colorbar
set(gca,'FontName','Latin Modern Roman','FontSize',22);
set(gca,'TickLabelInterpreter', 'latex');
axis square;ylabel('$y$ (m)','Interpreter','latex');xlabel('$x$ (m)','Interpreter','latex');
