function [rotation]=rodrigues(angle,axis)

axis=axis/norm(axis);
rotation=eye(3)+(sin(angle))*skew(axis)+(1-cos(angle))*(skew(axis))^2;


end