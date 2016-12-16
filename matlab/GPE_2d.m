%% Effective 2D single component Gross-Pitaevskii equation solver using split-operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear;
figure;
 %%Trapping frequencies in X,Y,Z dimensions. 2D is same as tightly confined
 %%in Z dimension.
 omegax = 2*pi*1;
 omegay=2*pi*1;
 omegaz = 2*pi*10; %Used for determining interaction strength of 2d system
 
 %%Constants required for simulation
 m=1.4431607e-25; %Rb87 mass
 global hbar; hbar=1.05457e-34;
 N=4.9e5;%Number of atoms
 as = 4.76e-9; %S-wave scattering length
 a0x = sqrt(hbar/(2*m*omegax))%Trapped width in X dim
 a0y = sqrt(hbar/(2*m*omegay))%Trapped width in Z dim
 
 Omega=0.25*omegax;%This defines the rotation frequency of trap.
 Rxy = (15.^0.2)*(N*4.76e-9*sqrt(m*omegaz/hbar)).^0.2;
 
 xmax=5*Rxy*a0x; ymax=5*Rxy*a0x; 
 Ngx=2^8; Ngy=2^8;
 [x,dx,px,dpx]=fftdef(xmax,Ngx);
 [y,dy,py,dpy]=fftdef(ymax,Ngy);
 
 U=(N)*(4*pi*hbar*hbar*as/m)*sqrt(m*(omegaz)/(2*pi*hbar)); %The interaction term I use. Yours will probably be different.
 
 %%%  grids in position/momentum space 
 [ym,xm]=meshgrid(y,x); [pym,pxm]=meshgrid(py,px);
 r=sqrt(xm.^2+ym.^2);   p=sqrt(pxm.^2+pym.^2); 
 xpyv=xm.*pym;          ypxv=pxm.*ym;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2pi phase singularity for vortex creation/annihilation
for iii=1:Ngx
for jjj=1:Ngy
if y(iii)>=0
PH(jjj,iii)=atan((y(jjj)+(dx/1000))./x(iii)) -pi/2;
elseif x(iii)<0
PH(jjj,iii)=atan((y(jjj)+(dx/1000))./x(iii)) +pi/2;
end
end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Phase winding
 l=0;
 %%%  starting wavefunction 
 wfc=exp(-( (xm./(Rxy*a0x)).^2 + (ym./(Rxy*a0y)).^2));%My initial gaussian wavefunction guess.
 wfc=wfc.*exp(1i*l*PH'); %Apply phase profile
 wfc=wfc/(sum(sum(abs(wfc).^2))*dx*dy)^(1/2); %Normalisation

 %%%  Hamiltonian operators
 
 V=1/2*m*(omegax^2.*xm.^2 + omegay^2.*ym.^2);
 K=hbar^2/(2*m)*(pxm.^2 + pym.^2);%Momentum
 
 % Optical lattice. Adjust spacing to see fit.
 clear A B C
 lattice_width=8*780e-9;
 x_shift=0;
 y_shift=0;
 [Xc,Yc]=meshgrid(x,y);
 A = cos((2*sqrt(3)*(xm+x_shift)/2 - (ym+y_shift)/2)/lattice_width).^2;
 B = cos((2*sqrt(3)*(xm+x_shift)/2 + (ym+y_shift)/2)/lattice_width).^2;
 C = cos(((ym+y_shift))/lattice_width).^2;
 
 A = cos((-(xm+x_shift)/2 + sqrt(3)*(ym+y_shift)/2)/lattice_width).^2;
 B = cos(((xm+x_shift)/2 + sqrt(3)*(ym+y_shift)/2)/lattice_width).^2;
 C = cos(((xm+x_shift))/lattice_width).^2;

 V_opt = A+ B+ C;
 figure;pcolor(V_opt);shading interp;axis square
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Omega=omegax*(ii/100);

 Dt=1e-4%Timestep for evolution
 GVr=(exp(-((  V + 0*(hbar*omegax).*V_opt')./hbar)*(Dt/2))); GVr2=exp(-(V./hbar)*Dt);  GKp=(exp(-(K./hbar)*Dt));  %Imaginary time unitary evolution operators for position, momentum
 GL2=exp(-Omega*(-ypxv)*Dt);  GL1=exp(-Omega*xpyv*Dt); %and angular momentum (for rotation only)

 UhbarDt=U./hbar*(Dt/2);%Just combining the interaction term with the form needed for time evolution
%% This performs imaginary time evolution to find the groundstate.
 for c1=1:500;
c1;
tic
for a=1:100
    a
   %Symmetric splitting takes error from dt^2 to dt^3 
   
   wfc=GVr.*exp(-abs(wfc).^2.*(UhbarDt)).*wfc;
   
   wfc=ifft2(GKp.*fft2(wfc));
   
   wfc=GVr.*exp(-abs(wfc).^2.*(UhbarDt)).*wfc;
   
   if l ~= 0
       if mod(a,2)==1
            wfc=ifft(GL1.*fft(wfc,[],2),[],2);     %%These two lines handle the angular momentum operator. You need to FFT in 1 dimension, evolve, then fft in the other dimension and the same. This is for 2D.
            wfc=ifft(GL2.*fft(wfc));
       else
            wfc=ifft(GL2.*fft(wfc));
            wfc=ifft(GL1.*fft(wfc,[],2),[],2);   %%These two lines handle the angular momentum operator. You need to FFT in 1 dimension, evolve, then fft in the other dimension and the same. This is for 2D.
       end
   end
   wfc=wfc./(sum(sum(abs(wfc).^2))*dx*dy)^(0.5); %Renormalise for groundstate

end
toc
%[E,L]=energy_ang(V,K,U,xm,ym,wfc,dx,dy);
%En(c1)=E; Lz(c1)=L;
pcolor(abs(wfc).^2);shading interp; axis square;
drawnow;
 end
 
 %% This performs real time evolution of the groundstate.
 
 Dt=1e-5%Timestep for real-time evolution
 EVr=exp(-i.*((  V + 0*(hbar*omegax).*V_opt')./hbar)*(Dt/2)); EVr2=exp(-i.*(V./hbar)*Dt);   EKp=exp(-i.*(K./hbar)*Dt);  %Imaginary time unitary evolution operators for position, momentum
 EL2=exp(-i.*Omega*(-ypxv)*Dt);  EL1=exp(-i.*Omega*xpyv*Dt); %and angular momentum (for rotation only)

 UhbarDt=U./hbar*(Dt/2);%Just combining the interaction term with the form needed for time evolution
 for c1=1:500;
    c1;

    for a=1:100
       %FFT to momentum space, and apply momentum evolution oeprator

       wfc=EVr.*exp(-abs(wfc).^2.*(UhbarDt)).*wfc;

       wfc=ifft2(EKp.*fft2(wfc));

       wfc=EVr.*exp(-abs(wfc).^2.*(UhbarDt)).*wfc;

       if mod(a,2)==1
            wfc=ifft(EL1.*fft(wfc,[],2),[],2);     %%These two lines handle the angular momentum operator. You need to FFT in 1 dimension, evolve, then fft in the other dimension and the same. This is for 2D.
            wfc=ifft(EL2.*fft(wfc));
       else
            wfc=ifft(EL2.*fft(wfc));
            wfc=ifft(EL1.*fft(wfc,[],2),[],2);   %%These two lines handle the angular momentum operator. You need to FFT in 1 dimension, evolve, then fft in the other dimension and the same. This is for 2D.
       end

    end

    %[E,L]=energy_ang(V,K,U,xm,ym,wfc,dx,dy);
    %En(c1+500)=E; Lz(c1+500)=L;
    pcolor(abs(wfc).^2);shading interp; axis square;
    drawnow;
end