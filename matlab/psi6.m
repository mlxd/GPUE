function [psi6_pos,nn] = psi6(pos,X,Y,radius)
%Calculate the orientational order parameter defined at point pos
%between the points (X,Y)
% pos:       Defines the location to calculate the orientational order
% X,Y:       Vector of entire X,Y range of points
% radius:    Radius over which to determine neighbouring points
%Returns
% psi6_pos:  The value of orientational order psi_6 at position pos
% nn:        Number of nearest neighbours

    [nn,idx] = findNN(pos,X,Y,radius); %find number of neighbours and indices of
                                        %neighbours for X,Y
    psi6_pos = 0;
    
    if size(nn,1) > 0
        for ii=1:size(nn,1)
            psi6_pos = psi6_pos + exp(6*1i*getAngle(pos, [X(idx(ii)), Y(idx(ii))] ));
        end
        psi6_pos = psi6_pos./length(nn);
    end
    nn=length(nn);
end
