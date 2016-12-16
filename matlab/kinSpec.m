function [EKc,EKi] = kinSpec(wfc_str,dimSize,dx,mass,qSpec,init,incr,fin)
%This function is used to process a large number of wavefunctions to examine
%the kinetic energy spectra. I'd love to port these to CUDA, but time is a
%precious resource that I have little of right now
%   wfc_str is the string of the wavefunction name. For time-ev this is
%       usually wfc_ev; wfc_0_const or wfc_0_ramp for ground-states
%   dimSize is a vector with the dimenion size, so that the wavefunctions
%       can be appropriately reshaped when loaded
%   dx is the increment along x. dx==dy
%   qSpec enables classical (0) or quantum (1) variant of energy spectra
%   init,incr,fin are the loop values for specifiying the range of the data

%Load initial wavefunction to setup grid
wfc = reshape(load([wfc_str,'_',int2str(0)]) + 1i*load([wfc_str,'i_',int2str(0)]),dimSize(1),dimSize(2));
[~,~,EKc,EKi,k_bin] = quKineticSpec(wfc,mass,dx,qSpec,1,0);
for ii=init:incr:fin
    wfc = reshape(load([wfc_str,'_',int2str(ii)]) + 1i*load([wfc_str,'i_',int2str(ii)]),dimSize(1),dimSize(2));
    [~,~,ekc,eki,~] = quKineticSpec(wfc,mass,dx,qSpec,1,ii);
    EKc=vertcat(EKc,ekc);EKi=vertcat(EKi,eki); %Prealloc may improve speed here
end
save Ek.mat;
