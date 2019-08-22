function plot_clusters(vGlyphs, T, bLines, bPlatonic)
% Show glyphs overlaid in clusters

	if(~exist('bLines'))
		bLines = false;
	end
	if(~exist('bPlatonic'))
		bPlatonic = true;
	end

	clf;
	nClusters = max(T);

	vColors = lapham( nClusters );
	for i = 1:nClusters
		subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
		hold on;
		iClusterGlyphs = find(T==i)';
		fAlpha = sqrt(1/length(iClusterGlyphs));
		for j = iClusterGlyphs
			if(bLines)
				plot(vGlyphs{j}(:,1), -vGlyphs{j}(:,2), '-', 'Color', vColors(i,:));
				plot(vGlyphs{j}(1,1), -vGlyphs{j}(1,2), '.', 'Color', vColors(i,:));
			else
				fill(vGlyphs{j}(:,1), -vGlyphs{j}(:,2), vColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', fAlpha);
			end
		end
		
		if(bPlatonic)
			[vMeanGlyph, vMedianGlyph] = platonicglyph(vGlyphs(iClusterGlyphs));
			plot(vMeanGlyph(:,1), -vMeanGlyph(:,2), '-', 'Color', [1 1 1]*0.99);
			plot(vMedianGlyph(:,1), -vMedianGlyph(:,2), '-', 'Color', [0 0 0]*0.99);
		end
		
		title(sprintf('%i',i));
		hold off;
		axis off;
		axis equal;
	end
	hold off;