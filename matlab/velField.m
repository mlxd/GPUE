function [] = velField(wfc0,dx,m,incr,x,y,normed, lims)
%velField calculates velocity field of wavefunction
%   velField plots the magnitude and direction of the velocity field of the
%   given wavefunction wfc0 in 2D.
%   dx is the increment along x, and assumes dx==dy
%   m is the mass of the atomic species
%   incr is the increment over which to calculate field direction. Too low
%       and it may be very dense, too high and very sparse. Start at 1, and
%       increase until happy
%   x,y are the grid spacings along the x and y axis
%   normed specifies whether to normalise the vector directions. 1 if yes,
%   0 otherwise
%   lims is [xMin xMax yMin yMax]. Hides the edge garbage. Otherwise, []

phase= angle(wfc0);
np1=unwrap(phase,[],1);
np2=unwrap(phase,[],2);

[velnp1x,velnp1y] = gradient(np1,dx,dx);
[velnp2x,velnp2y] = gradient(np2,dx,dx);

hbar=1.05457e-34;
v_y = (hbar/m)*(velnp1y);
v_x = (hbar/m)*(velnp2x);

valsx=[1:incr:size(wfc0,1)];%[200:1:300];
valsy=valsx;

%
if normed==1
	L=sqrt(v_x(valsy,valsx).^2 + v_y(valsy,valsx).^2);
else
    L=1;
end

imagesc(y(valsx),x(valsy),sqrt(v_x(valsy,valsx).^2+v_y(valsy,valsx).^2));hold on

q = quiver(y(valsx),x(valsy),v_x(valsy,valsx)./L,v_y(valsy,valsx)./L,'AutoScaleFactor',0.5,'Color','w');axis square;
hold off;

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',30);
xlabel('$x$ (m)');ylabel('$y$ (m)');
h=colorbar;set(h,'TickLabelInterpreter','latex')
ylabel(h,'|V| (m/s)')

if length(lims~=0)
    axis(lims);
end
