

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
load ../glyph_data01;

% Let's make extra sure that the glyphs are sorted before we get going
vGlyphs = sort_glyphs(vGlyphs);
v = cell2mat(vGlyphs(:,1:5));

% Cell array for table of data
vGlyphTable = cell(nGlyphs,6);
iTableIndex = 0;

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
			vGlyphsInWord = v(v(:,4)==iWord,5);				% DO NOT SORT
			
			% Plot and save image of glyph and add to table
			for l = 1:length(vGlyphsInWord)
				
				iGlyph = vGlyphsInWord(l);
				iGlyphIndex = find( v(:,5) == iGlyph );
				
				iTableIndex = iTableIndex+1;
				vGlyphTable(iTableIndex,:) = { iTableIndex, iGlyph, i, j, k, l };
				
				% Skip partial glyphs
				if(v(iGlyphIndex,1))
					vGlyphPoly = vGlyphs{iGlyphIndex,6};

					figure(l); clf; hold on;
					iLineWidth = 5;
					fill(vGlyphPoly(:,1), -vGlyphPoly(:,2), [1 1 1]/2, 'LineWidth', 5, 'EdgeColor', [0 0 0], 'FaceAlpha', 1);
					hold off;
					axis off;
					axis([ min(vGlyphPoly(:,1)) max(vGlyphPoly(:,1)) min(-vGlyphPoly(:,2)) max(-vGlyphPoly(:,2)) ] + ([-1 1 -1 1]*iLineWidth));
					axis equal;
					
					imgGlyph = print('-RGBImage');
					[h w d] = size(imgGlyph);
					
					imgGlyph = imresize(imgGlyph, 200.0/h);
					

					strImageFilename = sprintf('images%s%04i.png', filesep, iGlyph);
					strImageFilename = [ strWebFolder strImageFilename ];					% Replace this with a better system
					imwrite(imgGlyph, strImageFilename);
				end
			end
		end
	end
end


%%
tGlyphTable = cell2table(vGlyphTable);
writetable(tGlyphTable,'Glyph Table.csv');

% Save the data again, just in case anything gets changed here
save(['..' filesep 'glyph_data01b'],'vGlyphs','nGlyphs');




%% Push the website repo to upload the photos
% Replace this with a better system
% If you decide to use this system, make your own github.io page and change the path below

strCmd = 'cd /Users/christiancasey/Dropbox/Website/christiancasey.github.io; ';
strCmd = [ strCmd 'git add *; git commit -m "updating glyph images"; git push;' ];
[status,cmdout] = system(strCmd)



%% Script Termination
show_completion_message;

















































