

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap lapham;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');

iFig = 0;

%% Load working data from previous step
load glyph_data01;

%% Generate polygon vectors

iPtsHigh = 100;
iPtsLow = 20;

nTh = 24;
nR = 10;
maxR = log(get_glyph_extent(vGlyphs(:,6)));
mGlyphHistograms = zeros(nGlyphs,nTh*nR);

mGlyphPolygons = zeros(nGlyphs,iPtsLow*2);
vGlyphPolygons = cell(nGlyphs,1);

for i = 1:nGlyphs

	vGlyphPolygons{i,1} = vGlyphs{i,6};
	%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);

	% Interpolate the polygon so that all glyphs have the same # of points
	% Large interpolation (high point) - for display, etc.
	vGlyphPolygons{i,2} = interppoly(vGlyphPolygons{i,1}, iPtsHigh);
	% Small interpolation (low point) - for clustering
	vGlyphPolygons{i,3} = interppoly(vGlyphPolygons{i,1}, iPtsLow+1);
	vGlyphPolygons{i,3} = vGlyphPolygons{i,3}(1:end-1,:);

	% Subtract the centroid to zero out polygons and allow analysis
	vGPCenter = mean(vGlyphPolygons{i,1});
	for j = 1:2
		vGlyphPolygons{i,2}(:,j) = vGlyphPolygons{i,2}(:,j)-vGPCenter(j);
		vGlyphPolygons{i,3}(:,j) = vGlyphPolygons{i,3}(:,j)-vGPCenter(j);
	end
	vGlyphPolygons{i,4} = vGPCenter;
	
	% Add a centered polyshape for use in exact distance calculations
	warning('off','MATLAB:polyshape:repairedBySimplify');	% Temporarily disable warning, polyshape() should make corrections
	vGlyphPolygons{i,5} = polyshape(vGlyphPolygons{i,2}(:,1), vGlyphPolygons{i,2}(:,2), 'Simplify', true);
	warning('on','MATLAB:polyshape:repairedBySimplify');
	
	% Add a log-polar histogram to test in clustering
	m = logpolarhist(vGlyphs{i,6}, nTh, nR, maxR);
	vGlyphPolygons{i,6} = m;
	mGlyphHistograms(i,:) = m(:)';

	% Add the glyphs polygons to a big matrix
	mGlyphPolygons(i,:) = vGlyphPolygons{i,3}(:)';
end

%% Save working data for next step
save( 'glyph_data02', 'vGlyphs', 'nGlyphs', 'vGlyphPolygons', 'mGlyphPolygons', 'mGlyphHistograms' );

	
	
	
	
	
	
	
	
	
	
	
	








