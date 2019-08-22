function fMax = get_glyph_extent(v)
% Returns the maximum distance of a single point from the centroid in a set of glyphs.

	n = size(v,1);
	vMax = zeros(n,1);
	for i = 1:n
		v{i} = v{i} - mean(v{i});
		vMax(i) = max(sqrt(sum(v{i}.^2,2)));
	end
	fMax = max(vMax);