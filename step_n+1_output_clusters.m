

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');


%% Load working data from previous step
load glyph_data;

nGlyphs = size(vGlyphs,1);

%% Loop through clusters and format output
figure(401);
clf;


vPolyAverage = zeros(iPtsHigh,2);

for i = 1:nClusters
	iCluster = iSortedClusters(i);
	nGlyphsInCluster = nT(i);
	iClusterGlyphs = find(T==iCluster)';
	
	vClusterLineInstances = cell2mat( vGlyphs(iClusterGlyphs,3) );
	vClusterLines = unique(vClusterLineInstances);
	vNClustersPerLine = histc(vClusterLineInstances, vClusterLines);
	
	% Get average glyph per cluster
	for j = iClusterGlyphs
		vPolyAverage = vPolyAverage + vGlyphPolygons{j,2};
	end
	
	vPolyAverage = vPolyAverage / nGlyphsInCluster;
	
	subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
	plot(vPolyAverage(:,2),-vPolyAverage(:,1));
	axis equal
end
