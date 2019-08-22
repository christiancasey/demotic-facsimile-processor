function vR = getmagnitude(vXY)
% Returns the magnitudes from the centroid for a set of 2d vectors (x, y)

% 	% Center the vectors first
% 	vC = mean(vXY,1);
% 	if(~all(vC==0))
% 		for i = 1:2
% 			vXY(:,i) = vXY(:,i) - vC(i);
% 		end
% 	end
	
	vR = sqrt( sum( vXY.^2, 2 ) );
	