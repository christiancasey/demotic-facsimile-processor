function img = separate_connect(imgGlyphs, imgSeparators, imgConnectors)
% Returns a new image of glyphs with separations and connections included
	
	img = imgGlyphs;
	
	% Separation is easy, every place there is a separator, delete those pixels
	img(imgSeparators(:)) = 0;
	
	% Connection is more difficult.
	% If a connector does not connect any glyphs in this image, it should not be included, 
	% because it would be left floating in the image as a spurious tiny glyph.
	
	imgConnectors = bwlabeln(imgConnectors);
	nConnectors = max(imgConnectors(:));
	
	for i = 1:nConnectors
		
		imgI = (imgConnectors==i);
		
		if( any(img(imgI(:))) )
			img(imgI) = 1;
		end
	end
	