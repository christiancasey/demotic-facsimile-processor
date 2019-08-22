

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
whitebg('k');	colormap iris;
plottools('off');

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');

%%

clc

vTextImageParts = {
	'Page Breaks';
	'Line Breaks';
	'Word Breaks';
	'Morpheme Breaks';
	'Separators';
	'Connectors';
	'Full Glyphs 1';
	'Full Glyphs 2';
	'Full Glyphs 3';
	'Partial Glyphs 1';
	'Partial Glyphs 2';
	'Classifiers';
	'Lacunae';
	'Outline';
	'Facsimile';
	'Photograph';
};
dTextImageParts = containers.Map();

strFolder = 'input';
strExtension = 'png';

vImageSize = [];

for i = 1:length(vTextImageParts)
	
	strFilename = [strFolder filesep vTextImageParts{i} '.' strExtension];
	
	if isfile(strFilename)
		img = imread(strFilename);
		
		% For monochrome image parts, convert binary matrix
		if( i < min(find(strcmp(vTextImageParts(:,1),'Facsimile'))) )
			% All grayscale values map to 0, information must be some color
			img = get_mono_image(img);
			img = bwareaopen(img, 8);
		end
		
		% If there is a one pixel white border around the image, cut it out
		if( all( all( img(1,1,:) == img(2:end,1,:) ) ) )
			img = img(2:end-1,2:end-1,:);	% Cut out one-pixel border
		end
		
		% Make sure all images are identical in size
		vImageSizeI = size(img);
		if ~isempty(vImageSize)
			if( ~all( vImageSize(1:2) == vImageSizeI(1:2) ) )
				error('Images are not of equal size. Repeat layer export, adding a white one-pixel border if necessary.');
			end
		else
			vImageSize = vImageSizeI;
		end
		
		dTextImageParts(vTextImageParts{i}) = img;
	end
	
	% Display images
	if( isKey(dTextImageParts, vTextImageParts{i}) )
		figure(i);
		imagescz(img(:,:,1));
		title(vTextImageParts{i,1});
		getframe;
	end
end

vImageSize = vImageSize(1:2);

%% Pages, lines, sublines

% Mark page locations

if( ~isKey(dTextImageParts, 'Page Breaks') )
	dTextImageParts('Page Breaks') = ones(vImageSize);
end

imgPages = dTextImageParts('Page Breaks');

% If most of the image is zero,
% invert it so that blobs are ones and breaks are zeros
if( sum( imgPages(:) == 0 ) > sum( imgPages(:) == 1 ) )
	imgPages = 1-imgPages;
end

% Use conn comp and centroids to sort pages in order
% (Right-to-Left for Demotic, x values decreasing)
cc = bwconncomp(imgPages);
nPages = cc.NumObjects;
vPagesCentroid = regionprops(cc,'Centroid');
vPageCenters = zeros(nPages,2);
for j = 1:nPages
	vPageCenters(j,:) = vPagesCentroid(j).Centroid;
end
[~,iSort] = sort(vPageCenters(:,1),'Descend');
cc.PixelIdxList = cc.PixelIdxList(iSort);

imgPages = imgPages*0;
for j = 1:nPages
	imgPages(cc.PixelIdxList{j}) = j;
end

figure(101);
imagescz(imgPages);
datacursormode('on');

% Mark line locations
if( ~isKey(dTextImageParts, 'Line Breaks') )
	dTextImageParts('Line Breaks') = ones(vImageSize);
end

imgLines = dTextImageParts('Line Breaks');

% If most of the image is zero,
% invert it so that blobs are ones and breaks are zeros
if( sum( imgLines(:) == 0 ) > sum( imgLines(:) == 1 ) )
	imgLines = 1-imgLines;
end

% Add in page breaks
imgLines = imgLines & ( imgPages > 0 );

% Use conn comp and centroids to sort lines in order
% (top to bottom, y values increasing)
cc = bwconncomp(imgLines);
nLines = cc.NumObjects;
vLinesCentroid = regionprops(cc,'Centroid');
vLineCenters = zeros(nLines,2);
for j = 1:nLines
	vLineCenters(j,:) = vLinesCentroid(j).Centroid;
end
[~,iSort] = sort(vLineCenters(:,2),'Ascend');
cc.PixelIdxList = cc.PixelIdxList(iSort);


% Link each line to a page
vLinePages = zeros(cc.NumObjects, 1);

% Create labeled image
imgLines = imgLines*0;
for j = 1:nLines
	imgLines(cc.PixelIdxList{j}) = j;
	
	iLinePages = imgPages(cc.PixelIdxList{j});
	vUniquePages = unique( iLinePages );
	vNumPages = histc(iLinePages, vUniquePages);
	[~,iMax] = max(vNumPages);
	vLinePages(j) = vUniquePages(iMax);
end

figure(102);
imagescz(imgLines);
datacursormode('on');


% Mark subline locations
if( ~isKey(dTextImageParts, 'Word Breaks') )
	dTextImageParts('Word Breaks') = ones(vImageSize);
end

imgWords = dTextImageParts('Word Breaks');

% If most of the image is zero,
% invert it so that blobs are ones and breaks are zeros
if( sum( imgWords(:) == 0 ) > sum( imgWords(:) == 1 ) )
	imgWords = 1-imgWords;
end

% Add in page breaks
imgWords = imgWords & ( imgLines > 0 );

% Use conn comp and centroids to sort words in order
% (right to left, x values decreasing)
cc = bwconncomp(imgWords);
nWords = cc.NumObjects;
vWordsCentroid = regionprops(cc,'Centroid');
vWordCenters = zeros(nWords,2);
for j = 1:nWords
	vWordCenters(j,:) = vWordsCentroid(j).Centroid;
end
[~,iSort] = sort(vWordCenters(:,1),'Descend');
cc.PixelIdxList = cc.PixelIdxList(iSort);


% Link each subline to a line
vWordLines = zeros(cc.NumObjects, 1);

% Create labeled image
imgWords = imgWords*0;
for j = 1:nWords
	imgWords(cc.PixelIdxList{j}) = j;
	
	iWordLines = imgLines(cc.PixelIdxList{j});
	vUniqueLines = unique( iWordLines );
	vNumLines = histc(iWordLines, vUniqueLines);
	[~,iMax] = max(vNumLines);
	vWordLines(j) = vUniqueLines(iMax);
end

figure(103);
imagescz(imgWords);
datacursormode('on');



%% Process glyphs

disp('Creating glyph list...');

vGlyphKeys = { 'Full Glyphs 1', 'Full Glyphs 2', 'Full Glyphs 3', 'Partial Glyphs 1', 'Partial Glyphs 2' };
vFullGlyph = [ 1 1 1 0 0 ];
nGlyphSets = length(vGlyphKeys);
vImgGlyphs = cell(nGlyphSets,1);
vGlyphs = {};

% Get separator and connector layers
imgSeparators = zeros(vImageSize);
imgConnectors = zeros(vImageSize);
if( isKey(dTextImageParts, 'Separators') )
	imgSeparators = dTextImageParts('Separators');
end
if( isKey(dTextImageParts, 'Connectors') )
	imgConnectors = dTextImageParts('Connectors');
end

	

iGlyphCount = 0;
for i = 1:nGlyphSets
	
	if( ~isKey(dTextImageParts, vGlyphKeys{i}) )
		warning(sprintf('Image part does not exist: %s',vGlyphKeys{i}));
		continue;
	end
	disp(sprintf('Processing image part: %s', vGlyphKeys{i}));		% Progress status, because this is slow
		
	vImgGlyphs{i} = dTextImageParts(vGlyphKeys{i});
	
	% Not needed for classifier layer -- classifiers are already connected and separated
	vImgGlyphs{i} = separate_connect(vImgGlyphs{i}, imgSeparators, imgConnectors);
	
	vImgGlyphs{i} = bwmorph(vImgGlyphs{i}, 'spur');
	vImgGlyphs{i} = bwmorph(vImgGlyphs{i}, 'close', 4);
	
	% If most of the image is one,
	% invert it so that glyphs are ones and spaces are zeros
	if( sum( vImgGlyphs{i}(:) == 1 ) > sum( vImgGlyphs{i}(:) == 0 ) )
		vImgGlyphs{i} = 1-vImgGlyphs{i};
	end
	cc = bwconncomp(vImgGlyphs{i});
	nGlyphsI = cc.NumObjects;
	
	% REPLACED WITH sort_glyphs() below
% 	% Use centroids to sort glyphs (primarily for convenience)
% 	vGlyphsCentroid = regionprops(cc,'Centroid');
% 	vGlyphsCenters = zeros(nGlyphsI,2);
% 	for j = 1:nGlyphsI
% 		vGlyphsCenters(j,:) = vGlyphsCentroid(j).Centroid;
% 	end
% 	[~,iSort] = sort(vGlyphsCenters(:,2),'Descend');
% 	cc.PixelIdxList = cc.PixelIdxList(iSort);
% 	[~,iSort] = sort(vGlyphsCenters(:,1),'Ascend');
% 	cc.PixelIdxList = cc.PixelIdxList(iSort);
	
	% Create grid of pixel coordinates to use for glyph start point determination inside loop (find top left corner) 
	[h,w] = size(vImgGlyphs{i});
	
	% Create labeled image & link each glyph to a subline -> line -> page
	vImgGlyphs{i} = vImgGlyphs{i}*0;
	for j = 1:nGlyphsI
		iGlyphCount = iGlyphCount+1;
		
		vImgGlyphs{i}(cc.PixelIdxList{j}) = iGlyphCount;
		
		vXY = pixel_coordinates(cc.PixelIdxList{j}, h);
		[iStartX, iStartY] = topleftcorner(vXY);
		
		vGlyphPoly = bwtraceboundary(vImgGlyphs{i}==iGlyphCount, [iStartY, iStartX], 'N');
		vGlyphPoly = fliplr(vGlyphPoly);
		
		% Get the word to which each glyph belongs
		% Then also link it to the line and page
		iGlyphWords = imgWords(cc.PixelIdxList{j});
		vUniqueWords = unique( iGlyphWords );
		vNumWords = histc(iGlyphWords, vUniqueWords);
		[~,iMax] = max(vNumWords);
		iGlyphWord = vUniqueWords(iMax);
		iGlyphLine = vWordLines(iGlyphWord);
		iGlyphPage = vLinePages(iGlyphLine);
		
		vGlyphs(iGlyphCount,:) = { vFullGlyph(i), iGlyphPage, iGlyphLine, iGlyphWord, iGlyphCount, vGlyphPoly, vXY};
	end
end

%% Save working data for next step

% Sort glyphs by page, line, word, position in word (approximate, see sort_glyphs() documentation)
vGlyphs = sort_glyphs(vGlyphs);
nGlyphs = size(vGlyphs,1);

save('glyph_data01','vGlyphs','nGlyphs');



%% Script Termination
show_completion_message;


%% Go ahead and run table generation scripts

cd glyph_labels;
step_1b_create_glyph_table;
step_1c_create_word_table;








