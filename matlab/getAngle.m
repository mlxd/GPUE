function [ang] = getAngle(P1,P2)
% Returns angle between 2 points on a plane

ang = mod( atan2( P2(2) - P1(2),P2(1) - P1(1) ) , 2*pi );

end