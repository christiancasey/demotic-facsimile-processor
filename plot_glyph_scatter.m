function plot_glyph_scatter(X, vGlyphs, T)
% Show scatter plots of glyphs placed at point of first two dimensions of glyph vectors

	nClusters = max(T);
	nGlyphs = length(vGlyphs);
	
	
	[fMinDist, fMeanDist, fMedianDist, fMaxDist ] = get_distance_metrics(X);
	iScale = fMedianDist/get_glyph_extent(vGlyphs);
	
	
	clf; hold on;
	vColors = lapham(nClusters);
	for j = 1:nGlyphs
		fill(vGlyphs{j}(:,1)*iScale+X(j,1), -vGlyphs{j}(:,2)*iScale+X(j,2), vColors(T(j),:), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
	end
	axis equal;
	axis off;
	hold off;
	
	
	
	
	
	
	