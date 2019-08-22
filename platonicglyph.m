function [vMeanGlyph, vMedianGlyph] = platonicglyph(vGlyphs)
% Returns the platonic version of a set of overlapping glyphs

	nGlyphs = length(vGlyphs);
	if( ~iscell(vGlyphs) || (nGlyphs < 1) )
		error('Glyphs must be cell array of polygons');
	end
	
	nPoints = size(vGlyphs{1},1);
	
	mX = zeros(nGlyphs, nPoints);
	mY = zeros(nGlyphs, nPoints);
	
	for i = 1:nGlyphs
		
		vX = vGlyphs{i}(:,1);
		vY = vGlyphs{i}(:,2);
		
		mX(i,:) = vX(:)';
		mY(i,:) = vY(:)';
	end
	
	vMeanGlyph = [ mean(mX,1)' mean(mY,1)' ];
	vMedianGlyph = [ median(mX,1)' median(mY,1)' ];