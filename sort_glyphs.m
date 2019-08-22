function vGlyphs = sort_glyphs(vGlyphs)
% Returns a new version of vGlyphs with the glyphs in (approximately) the correct sequence
% NB At this time, the sorting does not account for stacked glyphs,
% because it is easy enough to correct these errors manually in the data entry step.
% This will be improved in a future version.

	% Sort by page, line, word, glyph
	v = cell2mat(vGlyphs(:,2:5));
	[v,i] = sortrows(v,1:4);
	vGlyphs = vGlyphs(i,:);
	v = cell2mat(vGlyphs(:,1:5));

	vGlyphsNew = cell(size(vGlyphs));
	iNewIndex = 0;

	% Dive down to the glyph level, keeping a running index for each part
	vPages = unique(v(:,2));
	for i = 1:length(vPages)

		iPage = vPages(i);
		vLinesInPage = unique(v(v(:,2)==iPage,3));
		for j = 1:length(vLinesInPage)

			iLine = vLinesInPage(j);
			vWordsInLine = unique(v(v(:,3)==iLine,4));
			for k = 1:length(vWordsInLine)

				iWord = vWordsInLine(k);
				vGlyphsInWord = unique(v(v(:,4)==iWord,5));
				nGlyphsInWord = length(vGlyphsInWord);

				% Use centroids to sort glyphs
				% Replace this with a sort that accounts for stacked glyphs
				vGlyphsCenters = zeros(nGlyphsInWord,2);
				for l = 1:nGlyphsInWord
					iGlyph = vGlyphsInWord(l);
					iGlyphIndex = find( v(:,5) == iGlyph );

					vGlyphsCenters(l,:) = mean(vGlyphs{iGlyphIndex,end}, 1);
				end
				[~,iSort] = sort(vGlyphsCenters(:,1),'Descend');
				vGlyphsInWord = vGlyphsInWord(iSort);

				% Plot and save image of glyph and add to table
				for l = 1:length(vGlyphsInWord)

					iGlyph = vGlyphsInWord(l);
					iGlyphIndex = find( v(:,5) == iGlyph );

					% Make a new vGlyphs array with glyphs in sort order
					iNewIndex = iNewIndex+1;
					vGlyphsNew(iNewIndex,:) = vGlyphs(iGlyphIndex,:);
				end
			end
		end
	end
	vGlyphs = vGlyphsNew;
	