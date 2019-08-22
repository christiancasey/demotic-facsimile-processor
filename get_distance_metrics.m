function [fMin, fMean, fMedian, fMax] = get_distance_metrics(x)
% Returns the distances between points in a set of vectors

	x = x(:,1:min(3,size(x,2)));
	n = size(x,1);
	k = size(x,2);
	
	% Old way finds all pairwise distances
% 	mDist = zeros(n);
% 	
% 	for i = 1:k
% 		mGrid = meshgrid(x(:,i));
% 		mDist = mDist + (mGrid-mGrid').^2;
% 	end
% 	
% 	mDist = sqrt(mDist);
% 	vDist = mDist(mDist~=0);

	% New way uses nearest neighbor
	iNN = knnsearch(x,x,'K',2);
	vDist = sqrt( sum( (x(iNN(:,2),:)-x).^2, 2 ) );
	
	fMin = min(vDist);
	fMean = mean(vDist);
	fMedian = median(vDist);
	fMax = max(vDist);