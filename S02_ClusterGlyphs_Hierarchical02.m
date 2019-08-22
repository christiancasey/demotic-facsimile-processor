

%% Script Initialization
switchToCD;
clc;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

%% Load image and begin processing

vFilenames = { 'Roman_Font_Test_Small.tif' };
vFilenames = { 'Roman_Font_Test.tif' };

strFilename = 'input/Classifiers.png';
% vFilenames = { 'Rosetta_Stone_All_Test.tif' };
% vFilenames = { 'Rosetta_Stone_All.tif' };
% vFilenames = { 'Rosetta_Stone_All.tif' 'Roman_Font_Test.tif' };


mGlyphPolygonsAll = [];

% for k = 1:length(vFilenames)

% 	strFilename = vFilenames{k}
	img = imread(strFilename);
	
	img = ( (img(:,:,1)==img(:,:,2)) & (img(:,:,2)==img(:,:,3)) );
	
% 	imshow(img);
	
% 	return;
	[iHeight iWidth] = size(img);

	% Convert to logical and invert so that bwlabel works properly
	img = ~(img > 0);

	% Remove very small blobs, since these are probably just noise
	img = bwareaopen(img, 500);
	%imagescz(img);

	% Flip the image horizontally and then transpose so that labels are in
	% approximately the right order. This is not a perfect solution.
% 	img = fliplr(img)';

	% Label blobs and get the number of distinct blobs for later looping
	%[imgGlyphs nGlyphs] = bwlabel(img);

	[vGlyphs,imgGlyphs] = bwboundaries(img,'noholes');
	nGlyphs = max(imgGlyphs(:));

	% Return both images to their orignal orientations.
% 	img = fliplr(img');
% 	imgGlyphs = fliplr(imgGlyphs');

	%% Display the labeled image and polygons
	%figure(1);
	%imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
	%hold on;

	iPtsHigh = 100;
	iPtsLow = 40;

	mGlyphPolygons = zeros(iPtsLow*2,nGlyphs);
	vGlyphPolygons = cell(nGlyphs,1);
	for i = 1:length(vGlyphs)
		
% 		% If the image was flipped, the polygon x values will be wrong
% 		vGlyphs{i}(:,1) = iWidth - vGlyphs{i}(:,1);
		
% 		vGlyphPolygon = vGlyphs{i};
		vGlyphPolygons{i,1} = vGlyphs{i};
		%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);

		% Interpolate the polygon so that all glyphs have the same # of points
		% Large interpolation (high point) - for display, etc.
		vGlyphPolygons{i,2} = interppoly(vGlyphPolygons{i,1}, iPtsHigh);
		% Small interpolation (low point) - for clustering
		vGlyphPolygons{i,3} = interppoly(vGlyphPolygons{i,1}, iPtsLow);


		% Use this plot function when the image is transposed
		%plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 1);
		% Use this plot to leave the original image unchanged
		%plot(vGlyphPolygon(:,2), vGlyphPolygon(:,1), 'w', 'LineWidth', 1);

		% Subtract the centroid to zero out polygons and allow analysis
		vGPCenter = mean(vGlyphPolygons{i,3});
		for j = 1:2
			vGlyphPolygons{i,3}(:,j) = vGlyphPolygons{i,3}(:,j)-vGPCenter(j);
		end
		vGlyphPolygons{i,4} = vGPCenter;

		%text(vGPCenter(1),vGPCenter(2),num2str(i),'Color','k','BackgroundColor','w');

		% Add the glyphs polygons to a big matrix
		mGlyphPolygons(:,i) = vGlyphPolygons{i,3}(:);
	end

	% Transpose mGlyphPolygons so that they're oriented properly 
	% [Glyphs Points]
	mGlyphPolygons = mGlyphPolygons';
	
	mGlyphPolygonsAll = [ mGlyphPolygonsAll ; mGlyphPolygons ];
% end

mGlyphPolygons = mGlyphPolygonsAll;


	
%% Get Principal Components and plot a subset of axes
mGPPC = mGlyphPolygons;
[mTPC,mGPPC,~] = getPC(mGlyphPolygons);
mTPC = eye(size(mTPC,1));
% mGPPC = mGlyphPolygons * mTPC;


%% Use some subset of PCs to create a clustering tree

iFig = 101;
close all;


X = mGPPC(:,1:20);
nClusters = 40

vDistMetrics = {
	'euclidean'					% 1
	'squaredeuclidean'			% 2
	'seuclidean'				% 3
	'cityblock'					% 4
	'minkowski'					% 5
	'chebychev'					% 6
	'mahalanobis'				% 7
	'cosine'					% 8
	'correlation'				% 9
	'spearman'					% 10
	'hamming'					% 11
	'jaccard'					% 12
};

vLinkMethods = {
	'single'	% 1		nearest distance (default)
	'complete'	% 2		furthest distance
	'average' 	% 3		unweighted average distance (UPGMA) (also known as group average)
	'weighted'	% 4		weighted average distance (WPGMA)
	'centroid'	% 5		unweighted center of mass distance (UPGMC)
	'median'	% 6		weighted center of mass distance (WPGMC)
	'ward'		% 7		inner squared distance (min variance algorithm)
};

for iDistMetric = [1 2 6 8 9]

	% X = mGlyphPolygons;
	Y = pdist(X,vDistMetrics{iDistMetric});
	sqY = squareform(Y);



	for iLinkMethod = [ 7 ]
		Z = linkage(Y, vLinkMethods{iLinkMethod});

		%figure(31);
		%dendrogram(Z,0);
		%zoom on;


		%%%%% Cluster
		T = cluster(Z,'maxclust',nClusters);

		% Sort clusters according to number of members and display the biggest
		% first
		nT = histc(T,1:nClusters);
		[nT,iSortedClusters] = sort(nT,'descend');

		% Show clustered glyphs in situ
		figure(iFig+100);
		iFig = iFig+1;
		clf;
		hold on;

		vColors = iris( nClusters );
		for i = iSortedClusters'

			iClusterGlyphs = find(T==i)';

			for j = iClusterGlyphs
				fill(vGlyphPolygons{j,2}(:,2), -vGlyphPolygons{j,2}(:,1), vColors(i,:));
				plot(vGlyphPolygons{j,3}(:,2)+vGlyphPolygons{j,4}(2), -vGlyphPolygons{j,3}(:,1)-vGlyphPolygons{j,4}(1), '-', 'Color', [1 1 1]/2);
				text(vGlyphPolygons{j,4}(2)+100,-vGlyphPolygons{j,4}(1)+100,sprintf('%d',i),'Color',vColors(i,:),'BackgroundColor',[1 1 1]*7/8)
			end

		end
		
		
		title(sprintf('Dist: %s, Link: %s (%i,%i)', vDistMetrics{iDistMetric}, vLinkMethods{iLinkMethod}, iDistMetric, iLinkMethod));
		axis equal;
		zoom on;
		getframe;
		hold off;
		zoom(4);

	end
end

%% Show glyphs on cluster map

return

clc
figure(201);
clf;
hold on;

iScale = 0.1;

vColors = iris( nClusters );
iStart = 1;%find(nT < 100, 1, 'first')
iEnd = nClusters%find(nT < 5, 1, 'first')
for i = iSortedClusters(iStart:iEnd)'
	
	iClusterGlyphs = find(T==i)';
	
	for j = iClusterGlyphs
		plot3d(vGlyphPolygons{j,3}(:,2)*iScale+mGPPC(j,1), -vGlyphPolygons{j,3}(:,1)*iScale+mGPPC(j,2), vGlyphPolygons{j,3}(:,1)*0 + mGPPC(j,3), '-', 'Color', vColors(i,:));
		
		text(mGPPC(j,1),mGPPC(j,2), mGPPC(j,3), sprintf('%d',i),'Color',vColors(i,:),'BackgroundColor',[1 1 1]*7/8)
	end
	
	getframe;
	%[ i n ]
end


