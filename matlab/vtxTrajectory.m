function [vorts] = vtxTrajectory(start, step, fin, initial,threeD, dx, dimSize, dt)
%vtxTrajectory trajectory plots of 2D vortex data.
%   vtxTrajectory takes the output of the postprocessed vortex positions
%   from the matlab script as loaded by loadVtx
%   Start step and fin are the initial, stepsize and final datasets
%   initial indicates if you wish to highlight the 
%   threeD also plots aspacetime (XYT) diagram in 3D.
%   dx is the stepsize along one dimension and assumes the same over x,y
%   dimSize is the length of the dimensions, and assumes x==y
%   dt is the time increment for the spacetime diagram.
% Note: Since the processed data assumes an initial starting position, the
% data at 0 is named vort_arr_0, while postprocessed data is vort_ord_xxx.
% Also, plotting assumes the data to be centred at 0 with for half dimSize
% of the system

    %Load vortices and convert struct to cell format
    vorts = struct2cell(loadVtx(start,step,fin,0));
    
    figure;
    for ii=1:size(vorts,2)
        if(ii==1) && (initial == 1)%load inital position and indicate with markers
            f0=csvread(['vort_ord_',int2str(start),'.csv'],1,0);
            p=scatter((f0(:,1)-dimSize/2)*dx,(f0(:,2)-dimSize/2)*dx,40,[0. 0. 0.],'filled','MarkerEdgeColor',[0. 0. 0.],'MarkerFaceColor',[0.4,0.4,0.4]); hold on;
            alpha(p,0.6);
        end
        plot(([vorts{1,ii,:}]-dimSize/2)*dx,([vorts{2,ii,:}]-dimSize/2)*dx,'color',rand(3,1),'LineWidth',2); 
        hold on;
    end

    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontName','Latin Modern Roman','FontSize',36);axis square;
    xlabel(' $x$ (m)','Interpreter','latex');ylabel( ' $y$ (m)','Interpreter','latex');
    print('-dpng','-r300',strcat('Trajectory.png',''));

    if threeD==1% Generate 2+1 spacetime plot of vortices
        figure
        for ii=1:length([vorts{1,:,1}])
            plot3(([vorts{1,ii,:}]-dimSize/2).*dx,([vorts{2,ii,:}]-dimSize/2).*dx,1e-2.*(1:length([vorts{2,ii,:}].*dt*step)));hold on
        end
        set(gca,'TickLabelInterpreter', 'latex');
        set(gca,'FontName','Latin Modern Roman','FontSize',22);
        xlabel('$x$ (m)','Interpreter','latex');ylabel('$y$ (m)','Interpreter','latex');zlabel('$t$ (s)','Interpreter','latex');
        print('-dpng','-r300',strcat('Trajectory3d.png',''));
    end
end