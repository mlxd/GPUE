function [neighbours,locations] = findNN(pos,X,Y,radius)
% Returns values and indices of nearest neighbours to location pos in
% positions X,Y within radius.

    r = sqrt(bsxfun(@minus,pos(1),X).^2 +  bsxfun(@minus,pos(2),Y).^2);
    locations = find(r < radius & r~=0);
    neighbours = [X(locations) Y(locations)];
end