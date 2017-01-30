function [psi6_pos,nn] = psi6_DT(pos,X,Y)
%Calculate the orientational order parameter defined at point pos
%between the points (X,Y)
% pos:       Defines the location to calculate the orientational order
% X:         Surrounding points in X coord to calculate oreintational order
% Y:         Surrounding points in Y coord to calculate oreintational order
%Returns
% psi6_pos:  The value of psi_6(pos)
% nn:        Number of nearest neighbours
    nn = length(X);
    psi6_pos = 0;
    
    if size(nn,1) > 0
        for ii=1:nn
            psi6_pos = psi6_pos + exp(6*1i*getAngle(pos, [X((ii)), Y((ii))] ));
        end
        psi6_pos = psi6_pos./nn;
    end
   
end
