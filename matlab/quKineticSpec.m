function [varargout] = quKineticSpec(Psir,m,dx,q,avg,iii)
% QUKINETICSPEC Evalutes the quantum kinetic energy spectrum [PRA 89, 053631 (2014)]
% Built on work by Dr A. White, and extended by (Dr?) Lee James O'Riordan.
%   Ek = quKineticSpec(Psir,m,dx,q,avg,iii), 
%   Psir N-dim wavefunction in r space with grid spacing dx and mass m
% 	q = include phase for quantum spectrum, 1=yes,0=no; 
%   avg = perform angled average. Likely always want this to be yes, 1
% 	iii output index for naming the figure files. Useful when looped with
% 	kinSpec script

% Setup
%##############################################################################%
hbar=1.05457e-34;
rDim = size(Psir); xDim = rDim(1); yDim = rDim(2);
%##############################################################################%
% Build camp in K-space
dkx = 2*pi/(dx*xDim); dky = 2*pi/(dx*yDim);
 
kx = [linspace(0,(xDim/2-1)*dkx,xDim/2) linspace(-xDim/2*dkx,-dkx,xDim/2)]';
ky = [linspace(0,(yDim/2-1)*dky,yDim/2) linspace(-yDim/2*dky,-dky,yDim/2)]';

[kxm,kym] = meshgrid(kx,ky);
%##############################################################################%
%|k|
km_mag = sqrt(kxm.^2 + kym.^2);
k_mag = sqrt(kx.^2 + ky.^2);
kMax = max(max(k_mag)); 

%##############################################################################%
Psik = fftn(Psir);

%##############################################################################%
%##############################################################################%
%   Velocity field calculation
%##############################################################################%

phase = angle(Psir);

np1=unwrap(phase,[],1);
np2=unwrap(phase,[],2);

[velnp1x,velnp1y] = gradient(np1,dx,dx);
[velnp2x,velnp2y] = gradient(np2,dx,dx);

v_y = (hbar/m)*(velnp1y);
v_x = (hbar/m)*(velnp2x);
%

v_x(find(v_x==0)) = 1e-100;
v_y(find(v_y==0)) = 1e-100;

u_x=exp(q*1i*angle(Psir)).*abs(Psir).*(v_x);
u_y=exp(q*1i*angle(Psir)).*abs(Psir).*(v_y);

%figure(4);pcolor(abs(u_x.^2 + u_y.^2));shading interp;title('|U|');colorbar;axis square

u_kx = fftn(u_x);
u_ky = fftn(u_y);

kxmkym = kxm.*kym;
uc_kx=(kxm.^2.*u_kx + kxm.*kym.*u_ky)./(km_mag.^2+1e-100);
uc_ky=(kxm.*kym.*u_kx + kym.^2.*u_ky)./(km_mag.^2+1e-100);

ui_kx = u_kx - uc_kx;
ui_ky = u_ky - uc_ky;

uc_x = ifft2(uc_kx); uc_y = ifft2(uc_ky);
ui_x = ifft2(ui_kx); ui_y = ifft2(ui_ky);

Ec = 0.5*abs(uc_x.^2 + uc_y.^2);
Ei = 0.5*abs(ui_x.^2 + ui_y.^2);

dimSize = size(Psir);
xx = linspace(-(dimSize(1)/2)*dx,(dimSize(1)/2)*dx,dimSize(1));
figure;pcolor(xx,xx,log10(Ec));shading interp;colorbar;axis square; title(['Comp ',int2str(iii)])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',24);
xlabel(' $x$ [m] ','Interpreter','latex')
ylabel(' $y$ [m] ','Interpreter','latex')

print('-dpng','-r300',['./Comp_CBAR_',int2str(iii),'.png']);
colorbar off;
print('-dpng','-r300',['./Comp_',int2str(iii),'.png']);

figure;pcolor(xx,xx,log10(Ei));shading interp;;colorbar;axis square;title(['Incomp ',int2str(iii)])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',24);
xlabel(' $x$ [m] ','Interpreter','latex')
ylabel(' $y$ [m] ','Interpreter','latex')

print('-dpng','-r300',['./Incomp_CBAR_',int2str(iii),'.png']);
colorbar off;
print('-dpng','-r300',['./Incomp_',int2str(iii),'.png']);

%figure;
[N1,C]=hist(abs(Psir.^2));drawnow;
hl = hbar/sqrt(m*2*(C(1))*1e6*(4*pi*hbar*hbar*4.67e-9/m)*sqrt(m*(100)/(2*pi*hbar)));

div=1;
ii=1;
for kk=1:length(k_mag)/2-1
    indX=find( k_mag(kk)<=km_mag & k_mag(kk+1) > km_mag);
    if avg==1
        div=length(indX);
    end
    ekc(ii) = (m*k_mag(kk)/(2)).*(sum(abs(sqrt(uc_kx(indX).^2 + uc_ky(indX).^2)).^2))./div;
    eki(ii) = (m*k_mag(kk)/(2)).*(sum(abs(sqrt(ui_kx(indX).^2 + ui_ky(indX).^2)).^2))./div;
    ii=ii+1;
end
figure;
loglog(k_mag(1:(xDim/2-1)),ekc,k_mag(1:(xDim/2-1)),eki,'LineWidth',2);%axis([1e4 1e7 5e-18 1e-10]);
%title('EK');%axis([1e4 1.2e7 1e-18 1e-6]);
legend({'$E^c(k)$','$E^i(k)$'},'FontSize',20,'FontWeight','bold','Interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'DefaultTextInterpreter','Latex')
set(gca,'FontName','Latin Modern Roman','FontSize',24);
xlabel(' $k$ [m$^{-1}$] ','Interpreter','latex')
ylabel(' $E(k)$ [a.u.] ','Interpreter','latex')
print('-dpng','-r300',['./EK_',int2str(iii),'.png']);
title(['Spectra ',int2str(iii)]);

varargout = {Ec,Ei,ekc,eki,k_mag};
%vline(2*pi/hl,'r','Healing');
%figure(4);
%loglog(k_mag(1:511),eki,'LineWidth',2);axis square;
%title('eki');


%ii=1;
%for kk=1:length(k_mag)/2-1
%    indX=find( k_mag(kk)<=km_mag & k_mag(kk+1) > km_mag);
%    Ek(ii) = (hbar^2*k_mag(kk)^3/(2*m))*sum(abs(Psik(indX)).^2)./length(indX);
%    ii=ii+1;
%end
%figure(5);
%loglog(k_mag(1:511),Ek);axis square;%axis([sqrt(2*dkx.^2) 2*pi/dx 1e-30 1e-5])
%title('Ek');
end
