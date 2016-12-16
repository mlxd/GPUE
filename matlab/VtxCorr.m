function [] = VtxCorr(dataSteps, dx, dimSize, radius, printIt, edgeCol)
%VtxCorr Plot the Voronoi diagram over the range of values spanned by
%dataSteps, and print the results to file. Essentially wraps
%voronoi2dCellColour and saves the resulting plots, with mean, variance and
%standard deviation values.
%
%dataSteps: The range of values to be plotted (e.g. [1e3 1e4 5e6])
%dx:        this is the spatial increment size of the data
%dimSize:   The number of elements along one dimension. Used to calculate
%               max and min spatial values. Assumes x==y
%radius:    defines the bounded region to calculate the Voronoi diagram and all 
%               quantities. Useful to avoid cells tending to infinity.
%               Assumes a circularly symmetric dataset.
%printIt:   1 if plots are to be saved, 0 otherwise.
%edgeCol:   Defines the Voronoi diagram edge colours for edge indexed pair. 
%               edgeCol==0 is black, edgeCol>0 is white
% Returns
%%Testcase: Generate diagrams for data at steps 1e3 1e4 and 1e5
% dx=1e-4; dimSize=1024; radius = 200*dx; printIt = 0; edgeCol = 1;
% VtxCorr([1e3 1e4 1e5], dx, dimSize, radius, printIt, edgeCol)

currentDirectory = pwd
[upperPath, deepestFolder, ~] = fileparts(currentDirectory) 

mArr=[];
sArr=[];
c=0
LB=flipud(lbmap(256,'RedBlue'));
for ii=[1e3 1e4 1e5 6e5]%
    hold on
    c=c+1;
    vtx=csvread(['vort_arr_',int2str(ii)],1,0);
    %vtx=csvread(['vort_ord_',int2str(ii),'.csv'],1,0); %Indexing needs to
    %be modified if you wish to use the ordered data sets. 
       
    %Calculate the Voronoi diagram of the resulting data, and plot with the
    %chosen color scheme. The local orientational correlations are returned
    %as p6cp6, and arr holds the values of the cell color scheme quantity.
    edges = edgeCol*ones(length(vtx(:,2)),1);
    [p6cp6, arr, ~,~,~]=voronoi2dCellColour((vtx(:,2)-(dimSize/2))*dx,(vtx(:,4)-(dimSize/2))*dx,radius,edges,dx,1)
    axis square;axis on;
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'DefaultTextInterpreter','Latex')
    set(gca,'FontName','Latin Modern Roman','FontSize',30);
    ylabel('$y$ (m)','Interpreter','latex');xlabel('$x$ (m)','Interpreter','latex');
    h=colorbar;set(h,'TickLabelInterpreter','latex');colormap (LB)
    
    p6cp6=p6cp6(find(p6cp6 ~= 0)); %Remove any 0 entries for sanity purposes
    
    %Determine mean, standard dev and variance values
    mArr(ii/1000 +1) = mean(p6cp6);
    vArr = var(arr);
    sArr(ii/1000+ 1) = std(p6cp6);
    if printIt %Naming scheme holds all useful info about data
        print('-depsc',['Voronoi_corr_',deepestFolder,'_t',int2str(ii),'_std',num2str(sArr(c)),'_m',num2str(mArr(c)),'.eps']);
    end
    drawnow;pause(0.5)
end


%% STD plot
clf
plot(linspace(0,6,length(sArr)),sArr,'LineWidth',2);set(gca,'TickLabelInterpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',30);
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$\sigma$','Interpreter','latex');
axis tight;set(gca,'PlotBoxAspectRatio',[1.0000    0.2613    0.2613]);
if printIt
    print('-depsc',['STD_',deepestFolder,'.eps']);
end
%%

%% Mean plot
clf
plot(linspace(0,6,length(mArr)),mArr,'LineWidth',2);set(gca,'TickLabelInterpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',30);
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$\bar{g_6}(0)$','Interpreter','latex');
axis tight;set(gca,'PlotBoxAspectRatio',[1.0000    0.2613    0.2613]);
if printIt
    print('-depsc',['AreaM_',deepestFolder,'.eps']);
end
