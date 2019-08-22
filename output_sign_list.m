function output_sign_list(strFilename, vGlyphLocations, vGlyphs, T)
% Create a sign list with platonic glyph, line numbers, and counts

	latex_wrapper_constants;		% Strings for LaTeX document
	
	strTeX = strTop;
	
	nClusters = max(T);

	vColors = lapham( nClusters );
	for i = 1:nClusters
		
		figure(i); clf; hold on;
		
		iClusterGlyphs = find(T==i)';
		nClusterGlyphs = length(iClusterGlyphs);
% 		fAlpha = sqrt(1/length(iClusterGlyphs));
		fAlpha = 0.6;
		for j = iClusterGlyphs
			fill(vGlyphs{j}(:,1), -vGlyphs{j}(:,2), vColors(i,:), 'EdgeColor', vColors(i,:), 'FaceAlpha', fAlpha);
		end
		
		[vMeanGlyph, vMedianGlyph] = platonicglyph(vGlyphs(iClusterGlyphs));
% 		plot(vMeanGlyph(:,1), -vMeanGlyph(:,2), '-', 'Color', [1 1 1]*0.99);
% 		plot(vMedianGlyph(:,1), -vMedianGlyph(:,2), '-', 'Color', [1 1 1]*0.99);
		fill(vMeanGlyph(:,1), -vMeanGlyph(:,2), [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5);
% 		fill(vMedianGlyph(:,1), -vMedianGlyph(:,2), [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5);
		
		hold off;
		axis off;
		axis equal;
		
		strImageFilename = sprintf('latex_output%sFigures%s%04i.eps', filesep, filesep, i);
% 		print(strImageFilename, '-depsc', '-painters');
		
		
		vLines = vGlyphLocations(iClusterGlyphs,3);
		
		vLinesUnique = unique(vLines);
		vLinesCounts = histc(vLines, vLinesUnique);
		vLinesOutput = [vLinesUnique(:) vLinesCounts(:)]';
		strLinesOutput = sprintf('%i (%i), ', vLinesOutput);
		if(length(strLinesOutput) < 2)
			error('Not enough instances of glyph for output');
		end
		
		strLinesOutput = strLinesOutput(1:end-2);
		strLinesOutput = strrep(strLinesOutput,' (1)',''); % Get rid of counts for 1
		
		strTableRow = sprintf('%i & \\includegraphics[align=t,width=0.8cm]{%04i} & %s & %i \\\\\n\\hline\n', ...
			i, i, strLinesOutput, nClusterGlyphs);
		strTeX = [ strTeX strTableRow ];
	end
	
	
	strTeX = [ strTeX strBtm ];
	txtwrite( strTeX, [ 'latex_output' filesep strFilename ] );
	
	
	