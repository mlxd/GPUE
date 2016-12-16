function [g6C] = g6_struct(X,Y,rad)
%Determine the orientational correlation function between every pairing of
%points, and sorts the results based on distance between points. This gives
%g6(r).
% X,Y:   Vectors of x,y values for points.
% rad:   Radius in which to examine for neighbouring points.
%Return:
% g6C:   Matrix of g6 values and parameters

%Determine all orientational order values for points
psi6p = zeros(length(X),1);
for ii=1:length(X)
    psi6p(ii) = psi6([X(ii),Y(ii)],X,Y,rad);
end

%Calculate all n choose k pairings of points, and give distance between
%them, values for the respective orientational orders, and g6 values. Sort
%the values based on separated distance
S = uniqPairIdx_precalc([X,Y],psi6p);
[~,order] = sort([S(:).rabs],'ascend');
SSorted = S(order);

%Create output structure for next stage of calculation
for ii=1:size(SSorted,2)
    g6S(ii).r = SSorted(ii).rabs;
    psi6_1 = SSorted(ii).cor0;
    psi6_2 = SSorted(ii).cor1;
    g6S(ii).val = abs(conj(psi6_1)*conj(psi6_1));
end
g6C = cell2mat(squeeze(struct2cell(g6S)));
end