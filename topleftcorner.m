function [x, y] = topleftcorner(vXY)
% Get the top left corner of a blob from a set of vectors (x,y)

	vC = mean(vXY,1);

	[~,iMin] = min( sum( (vC-vXY-1000).^2, 2 ) );
	x = vXY(iMin,1);
	y = vXY(iMin,2);