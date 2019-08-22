function vTh = getangle(vXY)
% Returns the angles from the centroid for a set of 2d vectors (x, y)

% 	% Center the vectors first
% 	vC = mean(vXY,1);
% 	if(~all(vC==0))
% 		for i = 1:2
% 			vXY(:,i) = vXY(:,i) - vC(i);
% 		end
% 	end
	
	vTh = atan(vXY(:,2)./vXY(:,1));
	vTh(vTh<0) = vTh(vTh<0)+pi;				% Rotate negative angles to positive
	vTh(vXY(:,2)<0) = vTh(vXY(:,2)<0) + pi;		% Rotate negative y values around
	