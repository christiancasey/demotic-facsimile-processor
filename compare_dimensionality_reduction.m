

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');



%% Compare manual counts

vCounts = [ 45 162 73 17 12 81 3 6 1 42 18 28 13 9 8 6 27 3 21 2 3 3 10 1 10 11 2 1 1 1 4 4 3 1 5 3 2 6 1 1 1 1 1 1 1 1 1 3 1 1 1 10 1 1 2 1 3 1 ];
vCounts = vCounts';
vCounts = sort(vCounts, 'descend');
sum(vCounts)


%% Load image and begin processing


strFilename = 'input/Classifiers.png';



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
	iPtsLow = 20;

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
		vGlyphPolygons{i,3} = interppoly(vGlyphPolygons{i,1}, iPtsLow+1);
		vGlyphPolygons{i,3} = vGlyphPolygons{i,3}(1:end-1,:);

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
		
		vGlyphPolygons{i,5} = polyshape(vGlyphPolygons{i,2}(:,2)-vGPCenter(2), vGlyphPolygons{i,2}(:,1)-vGPCenter(1));

		%text(vGPCenter(1),vGPCenter(2),num2str(i),'Color','k','BackgroundColor','w');

		% Add the glyphs polygons to a big matrix
		mGlyphPolygons(:,i) = vGlyphPolygons{i,3}(:);
	end

	% Transpose mGlyphPolygons so that they're oriented properly 
	% [Glyphs Points]
	tic
	mGlyphPolygons = mGlyphPolygons';
	toc
	mDist = polydist( vGlyphPolygons(:,5) );
	toc

	return;
	
%% Get average glyph width
vWidths = zeros(nGlyphs,1);
for i = 1:nGlyphs
	vWidths(i) = max(vGlyphPolygons{i,1}(:,1))-min(vGlyphPolygons{i,1}(:,1));
end
fGlyphWidth = median(vWidths);

%% Dimensionality Reduction

vDRMethods = {
	'PCA';				% 01	none
	'LDA';				% 02	none
	'MDS';				% 03	none
	'ProbPCA';			% 04	<int> max_iterations -> default = 200
	'FactorAnalysis';	% 05	none
	'GPLVM';			% 06	<double> sigma -> default = 1.0
	'Sammon';			% 07	none
	'Isomap';			% 08	<int> k -> default = 12
	'LandmarkIsomap';	% 09	<int> k -> default = 12
						%		<double> percentage -> default = 0.2
	'LLE';				% 10	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'Laplacian';		% 11	<int> k -> default = 12
						%		<double> sigma -> default = 1.0
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'HessianLLE';		% 12	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'LTSA';				% 13	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'MVU';				% 14	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'CCA';				% 15	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'LandmarkMVU';		% 16	<int> k -> default = 5
	'FastMVU';			% 17	<int> k -> default = 5
						%		<logical> finetune -> default = true
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'DiffusionMaps';	% 18	<double> t -> default = 1.0
						%		<double> sigma -> default = 1.0
	'KernelPCA';		% 19	<char[]> kernel -> {'linear', 'poly', ['gauss']} 
						%		kernel parameters: type HELP GRAM for info
	'GDA';				% 20	<char[]> kernel -> {'linear', 'poly', ['gauss']} 
						%		kernel parameters: type HELP GRAM for info
	'SNE';				% 21	<double> perplexity -> default = 30
	'SymSNE';			% 22	<double> perplexity -> default = 30
	'tSNE';				% 23	<int> initial_dims -> default = 30
						%		<double> perplexity -> default = 30
	'LPP';				% 24	<int> k -> default = 12
						%		<double> sigma -> default = 1.0
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'NPE';				% 25	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'LLTSA';			% 26	<int> k -> default = 12
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'SPE';				% 27	<char[]> type -> {['Global'], 'Local'}
						%		if 'Local': <int> k -> default = 12
	'Autoencoder';		% 28	<double> lambda -> default = 0
	'LLC';				% 29	<int> k -> default = 12
						%		<int> no_analyzers -> default = 20
						%		<int> max_iterations -> default = 200
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'ManifoldChart';	% 30	<int> no_analyzers -> default = 40
						%		<int> max_iterations -> default = 200
						%		<char[]> eig_impl -> {['Matlab'], 'JDQR'}
	'CFA';				% 31	<int> no_analyzers -> default = 2
						%		<int> max_iterations -> default = 200
	'NCA';				% 32	<double> lambda -> default = 0.0
	'MCML';				% 33	none
	'LMNN'				% 34	<int> k -> default = 3
};

nDim = 10;

vDRX = cell(34,2);

vBroken = [];
for iDRMethod = 1%:34%[1:11 18:29 31:34]
	
	disp(vDRMethods{iDRMethod});
	
	X = mGlyphPolygons;
% 	X = prewhiten(X);
	
	try
    	[X, mapping] = compute_mapping(X, vDRMethods{iDRMethod}, nDim);
		vDRX{iDRMethod,1} = X;
		vDRX{iDRMethod,2} = mapping;
	catch
		warning(sprintf('Error with method %i - %s', iDRMethod, vDRMethods{iDRMethod}));
		vBroken = [ vBroken iDRMethod ];
		continue;
	end
end


%% Clustering


nClusters = 64;

iFig = 100;
close all;

for iDRMethod = 1%[1:11 18:29 31:34]
	
% 	if any(vBroken==iDRMethod)
% 		continue;
% 	end
	X = mGlyphPolygons;
% 	X = abs(vDRX{iDRMethod,1});
% 	
% 	
% 	% Find minimum x distance in order to scale displayed glyphs
% 	[m1 m2] = meshgrid(X(:,1));
% 	mDist = abs(m1-m2);
% 	fSpace = median(mDist(mDist>0));
% 	iScale = fSpace/fGlyphWidth * 0.1;
	iScale = 0.1;
	
	% Use some subset of PCs to create a clustering tree

	vDistMetrics = {
		'euclidean'				% 1
		'squaredeuclidean'		% 2
		'seuclidean'			% 3
		'cityblock'				% 4
		'minkowski'				% 5
		'chebychev'				% 6
		'mahalanobis'			% 7
		'cosine'				% 8
		'correlation'			% 9
		'spearman'				% 10
		'hamming'				% 11
		'jaccard'				% 12
	};

	vLinkMethods = {
		'single'	% 1	nearest distance (default)
		'complete'	% 2	furthest distance
		'average' 	% 3	unweighted average distance (UPGMA) (also known as group average)
		'weighted'	% 4	weighted average distance (WPGMA)
		'centroid'	% 5	unweighted center of mass distance (UPGMC)
		'median'	% 6	weighted center of mass distance (WPGMC)
		'ward'		% 7	inner squared distance (min variance algorithm)
	};


% 	(1,2) (1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (3,4) (3,5) (4,5)
	
	
	% Former method combos: (1,4), (1,3), (6,3), (8,5)
	for iDistMetric = 1%[ 1 2 4 5 8 9 ]

% 		disp(vDRMethods{iDRMethod});
		
% 		Y = pdist(X,vDistMetrics{iDistMetric});
		Y = mDist;

		for iLinkMethod = 7%[ 3 4 5 6 ]

			iFig = iFig+1;


			%%%%% Cluster
			Z = linkage(Y, vLinkMethods{iLinkMethod});
			T = cluster(Z,'maxclust',nClusters);

			% Sort clusters according to number of members and display the biggest first
			nT = histc(T,1:nClusters);
			[nT,iSortedClusters] = sort(nT,'descend');

			% Show clustered glyphs in situ
			figure(iFig);
			clf; hold on;
			
			vColors = iris( nClusters );
			for i = iSortedClusters'
	
				iClusterGlyphs = find(T==i)';
	
				for j = iClusterGlyphs
					fill(vGlyphPolygons{j,2}(:,2), -vGlyphPolygons{j,2}(:,1), vColors(i,:));
					plot(vGlyphPolygons{j,3}(:,2)+vGlyphPolygons{j,4}(2), -vGlyphPolygons{j,3}(:,1)-vGlyphPolygons{j,4}(1), '-', 'Color', [1 1 1]/2);
% 					text(vGlyphPolygons{j,4}(2),-vGlyphPolygons{j,4}(1), 1000, sprintf('%d',i),'Color',vColors(i,:),'BackgroundColor',[1 1 1]*7/8)
				end
	
			end
			title(sprintf('Method No. %02i: %s', iDRMethod, vDRMethods{iDRMethod}));
			axis equal;
			zoom on;
			hold off;
			getframe;
% 			zoom(4);
			iEnd = nClusters;

% 			% Show clustered glyphs with first three dimensions
% 			figure(iFig+100);
% 			clf;
% 			hold on;
% 
% 			iEnd = nClusters;%5
% 			vColors = iris( iEnd );
% 			for i = 1:iEnd
% 				iCluster = iSortedClusters(i);
% 				iClusterGlyphs = find(T==iCluster)';
% 				for j = iClusterGlyphs
% % 					plot3d(vGlyphPolygons{j,3}(:,2)*iScale+X(j,1), -vGlyphPolygons{j,3}(:,1)*iScale+X(j,2), vGlyphPolygons{j,3}(:,1)*0 + X(j,3), '-', 'Color', vColors(i,:));
% 					fill3d(vGlyphPolygons{j,3}(:,2)*iScale+X(j,1), -vGlyphPolygons{j,3}(:,1)*iScale+X(j,2), vGlyphPolygons{j,3}(:,1)*0 + X(j,3), vColors(i,:));
% 					% 				text(X(j,1),X(j,2),X(j,3)+100, sprintf('%d',i),'Color',vColors(i,:),'BackgroundColor',[1 1 1]*7/8)
% 				end
% 				title(sprintf('Method No. %02i: %s', iDRMethod, vDRMethods{iDRMethod}));
% % 				axis equal;
% 				getframe;
% 			end
% 			hold off;
% 			getframe;
			
			% Show glyphs overlaid in clusters
			figure(iFig+200);
			clf;

			vColors = iris( iEnd );
			for i = 1:iEnd
				subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
				hold on;
				iCluster = iSortedClusters(i);
				iClusterGlyphs = find(T==iCluster)';
				for j = iClusterGlyphs
					plot(vGlyphPolygons{j,2}(:,2)-vGlyphPolygons{j,4}(2), -vGlyphPolygons{j,2}(:,1)+vGlyphPolygons{j,4}(1), '-', 'Color', vColors(i,:));
				end
				hold off;
				title(sprintf('Method No. %02i: %s', iDRMethod, vDRMethods{iDRMethod}));
% 				axis equal;
				getframe;
			end
			hold off;

			getframe;

		end
	end
	
end

disp(vBroken);











