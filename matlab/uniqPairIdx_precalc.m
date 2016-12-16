function [S] = uniqPairIdx_precalc(R,psi6p)
%Calculate all unique pairings of the orientational correlations, and
%create a struct to return that holds the distance between paired elements,
%and the respective orientational order values.
% R:        Vectors of R=[X Y] values for points.
% psi6p:    Orientational order values defined over the range of points
%Return:
% S:        Struct of orientational corder values cor0,cor1 and distance
%            between elements
    count = 0;
    for ii=1:(size(R,1)-1)
        for jj=ii+1:size(R,1)
            count = count +1;
            S(count).cor0 = psi6p(ii);
            S(count).cor1 = psi6p(jj);
            S(count).rabs = sqrt( sum( (R(ii,:) - R(jj,:) ).^2 ) );
        end
    end
end