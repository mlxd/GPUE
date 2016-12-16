function [psi6_pos,nn] = psi6(pos,X,Y,radius)
    [nn,idx] = findNN(pos,X,Y,radius);
    psi6_pos = 0;
    
    if size(nn,1) > 0
        for ii=1:size(nn,1)
            psi6_pos = psi6_pos + exp(6*1i*getAngle2(pos, [X(idx(ii)), Y(idx(ii))] ));
        end
        psi6_pos = psi6_pos./length(nn);
    end
    nn=length(nn);
end
