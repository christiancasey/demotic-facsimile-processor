

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap lapham;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');

iFig = 0;

%% Folder info 
% This is for saving images and uploading them to public server.
% The simple data entry solution (hack) used in this project is to put all of the Demotic signs
% in a Google sheet (because it supports images in cells)
% and then add the relevant info in adjacent cells so that it can be used later.
% First the images are saved in a folder that is part of my website,
% then the website is updated on github.io,
% then the images are obtained from predictable URLs by the Google sheet: "=IMAGE(URL)".
% There are many ways in which this might be replaced by a more elegant solution,
% but for now this gets the job done quickly and relatively effeciently.
% If you decide to use this system, create a github.io website from a folder on your local machine
% and change the folder paths in this file to your own.

strWebFolder = '/Users/christiancasey/Dropbox/Website/christiancasey.github.io/glyph_labels/';

%% Load working data from previous step
load ../glyph_data01b;

% Let's make extra sure that the glyphs are sorted before we get going
vGlyphs = sort_glyphs(vGlyphs);
v = cell2mat(vGlyphs(:,1:5));

% Cell array for table of data
vGlyphTable = cell(nGlyphs,5);
iTableIndex = 0;

% Dive down to the word level, keeping a running index for each part
vPages = unique(v(:,2));
for i = 1:length(vPages)
	
	iPage = vPages(i);
	vLinesInPage = unique(v(v(:,2)==iPage,3));
	for j = 1:length(vLinesInPage)
		
		iLine = vLinesInPage(j);
		vWordsInLine = unique(v(v(:,3)==iLine,4));
		for k = 1:length(vWordsInLine)
			
			iWord = vWordsInLine(k);
			vGlyphsInWord = v(v(:,4)==iWord,5);				% DO NOT SORT
			

			strImageFilename = sprintf('word_images%s%04i.png', filesep, iWord);
			strImageFilename = [ strWebFolder strImageFilename ];
			
			figure(k); clf; hold on;
			iLineWidth = 5;
			xMin = 1e10; yMin = 1e10;
			xMax = -1e10; yMax = -1e10;
			
			% Plot each glyph
			for l = 1:length(vGlyphsInWord)
				
				iGlyph = vGlyphsInWord(l);
				iGlyphIndex = find( v(:,5) == iGlyph );

				vGlyphPoly = vGlyphs{iGlyphIndex,6};

				fill(vGlyphPoly(:,1), -vGlyphPoly(:,2), [1 1 1]/2, 'LineWidth', 5, 'EdgeColor', [0 0 0], 'FaceAlpha', 1);
				
				xMin = min( xMin, min(vGlyphPoly(:,1)) );
				xMax = max( xMax, max(vGlyphPoly(:,1)) );
				yMin = min( yMin, min(-vGlyphPoly(:,2)) );
				yMax = max( yMax, max(-vGlyphPoly(:,2)) );
			end

			hold off;
			axis off;
			axis([ xMin xMax yMin yMax ] + ([-1 1 -1 1]*iLineWidth));
			axis equal;
			
			% Save image of word and add to table
			imgGlyph = print('-RGBImage');
			[h w d] = size(imgGlyph);

			imgGlyph = imresize(imgGlyph, 200.0/h);
			imwrite(imgGlyph, strImageFilename);

			iTableIndex = iTableIndex+1;
			vGlyphTable(iTableIndex,:) = { iTableIndex, iWord, i, j, k };
		end
	end
end


%%

tGlyphTable = cell2table(vGlyphTable);
writetable(tGlyphTable,'Word Table.csv');

% Save the data again, just in case anything gets changed here
save(['..' filesep 'glyph_data01c'],'vGlyphs','nGlyphs');

%% Push the website repo to upload the photos
% Replace this with a better system
% If you decide to use this system, make your own github.io page and change the path below

strCmd = 'cd /Users/christiancasey/Dropbox/Website/christiancasey.github.io; ';
strCmd = [ strCmd 'git add *; git commit -m "updating glyph images"; git push;' ];
[status,cmdout] = system(strCmd)



%% Script Termination
show_completion_message;

















































