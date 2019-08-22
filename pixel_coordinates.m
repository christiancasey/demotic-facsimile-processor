function vXY = pixel_coordinates(vP, height)
% Returns a 2d vector of pixel x and y coordinates based on a pixel list and the height of the image
	vP = vP(:);
	
	vXY = [ ceil(vP/height) mod(vP-1,height)+1 ];