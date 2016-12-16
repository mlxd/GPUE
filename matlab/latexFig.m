function latexFig(gca, fontSize, cbar, xtick, ytick, xticklabels, yticklabels)
%latexFig LaTeXify the plot axes
%   This sets the plotting axis on supplied figure to LaTeX formatted
%   {x,y}tick can be used to overwrite the automatically defined tick
%   locations
%   {x,y}ticklabels overwrites the labels at the positions defined by
%   {x,y}tick values
%   This assumes that the Latin Modern fonts are installed
%   Fontsize sets the size of the text
%   cbar enables the colorbar if needed

%   example:
%   pcolor(rand(10,10));
%   latexFig(gca,30,1,2:2:10,1:2:10,{'a1','b1','c1','d1','e1'},{'a2','b2','c2','d2','e2'})

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex');
            
if nargin >= 2        
    set(gca,'FontName','Latin Modern Roman','FontSize',fontSize);
end

if cbar == 1
    h=colorbar; set(h,'TickLabelInterpreter','latex');
end

if nargin >= 5
    set(gca,'XTick',xtick);
    set(gca,'YTick',ytick);
end

if nargin >= 7
    set(gca,'XTickLabels',xticklabels);
    set(gca,'YTickLabels',yticklabels);
end

%set(gcf,'PaperUnits','normalized');
%set(gcf,'PaperPosition',[0 0 1 1]);
%set(gca,'YTickLabels',{'\fontsize{14}{0}0','\fontsize{14}{0}$\frac{\pi}{24}$','\fontsize{14}{0}$\frac{\pi}{12}$','\fontsize{14}{0}$\frac{\pi}{8}$','\fontsize{14}{0}$\frac{\pi}{6}$'})

end