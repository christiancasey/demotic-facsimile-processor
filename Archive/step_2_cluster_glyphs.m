

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap lapham;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');



%% Compare manual counts

vCounts = [ 45 162 73 17 12 81 3 6 1 42 18 28 13 9 8 6 27 3 21 2 3 3 10 1 10 11 2 1 1 1 4 4 3 1 5 3 2 6 1 1 1 1 1 1 1 1 1 3 1 1 1 10 1 1 2 1 3 1 ];
vCounts = vCounts';
vCounts = sort(vCounts, 'descend');
sum(vCounts)


%% Load working data from previous step

load glyph_data;

nGlyphs = size(vGlyphs,1);

	% Return both images to their orignal orientations.
% 	img = fliplr(img');
% 	imgGlyphs = fliplr(imgGlyphs');


%% Get glyph maximum extents
vMax = zeros(nGlyphs,2);
for i = 1:nGlyphs
	for j = 1:2
		vMax(i,j) = max(vGlyphs{i,5}(:,j))-min(vGlyphs{i,5}(:,j));
	end
end
vMax = max(vMax);
iMaxDim = max(vMax);

%% Generate polygon vectors

iPtsHigh = 100;
iPtsLow = 20;

nTh = 24;
nR = 10;
maxR = log(iMaxDim/2+10);
mGlyphHistograms = zeros(nGlyphs,nTh*nR);

mGlyphPolygons = zeros(nGlyphs,iPtsLow*2);
vGlyphPolygons = cell(nGlyphs,1);

for i = 1:nGlyphs

	vGlyphPolygons{i,1} = vGlyphs{i,5};
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
	vGlyphPolygons{i,5} = polyshape(vGlyphPolygons{i,2}(:,1), vGlyphPolygons{i,2}(:,2));
	
	
	% Add a log-polar histogram to test in clustering
	m = logpolarhist(vGlyphs{i,6}, nTh, nR, maxR);
	vGlyphPolygons{i,6} = m;
	mGlyphHistograms(i,:) = m(:)';

	% Add the glyphs polygons to a big matrix
	mGlyphPolygons(i,:) = vGlyphPolygons{i,3}(:)';
end

%% Brute force dissimilarity measure (using polygons directly)
% Transpose mGlyphPolygons so that they're oriented properly 
% [Glyphs Points]
tic
% mGlyphPolygons = mGlyphPolygons';
toc
mDist = polydist( vGlyphPolygons(:,5) );
toc

%% Try clustering with different types of vectors to see which one works best

for iVectorType = 1:2

	%% Evaluate clustering algorithms
	nClusters = 64;
	
	if(iVectorType==1)
		X = mGlyphPolygons;	strTitle = 'Cluster Evaluation Using Polygons'; strVectorType = 'Polygon';
	elseif(iVectorType==2)
		X = mGlyphHistograms;	strTitle = 'Cluster Evaluation Using Log-Polar Histogram'; strVectorType = 'Histogram';
	end

	if ~(iVectorType== 3)

		vNumClusters = 1:64;
		vClusterEval = evalclusters(X,'linkage','Silhouette','KList', vNumClusters);

		[fEvalMax,iEvalMax] = max(vClusterEval.CriterionValues);
		nBestCluster = vNumClusters(iEvalMax)

		figure(21);
		clf; hold on;
		plot(vNumClusters, vClusterEval.CriterionValues);
		plot(nBestCluster, fEvalMax, '*r');
		title(strTitle);
		hold off;

		% print(imagefilename('Latin_Cluster_Eval'), '-depsc', '-painters');
		savewipimage(sprintf('Demotic_Classifier_Cluster_Eval (Vector Type - %s)', strVectorType));
	end
	
	
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

	for iDRMethod = 1:34

		if any(vBroken==iDRMethod)
			continue;
		end
		X = abs(vDRX{iDRMethod,1});

		iFig = 100;
		close all;
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

		% Former method combos: (1,4), (1,3), (6,3), (8,5)
		for iDistMetric = 2%[ 1 2 4 5 8 9 ]

			if(any(iVectorType==[1 2]))
				Y = pdist(X,vDistMetrics{iDistMetric});
			elseif(iVectorType==3)
				Y = mDist;	strVectorType = 'Brute Force'			% Brute force distance matrix
			end

			for iLinkMethod = 1:7

				iFig = iFig+1;


				%%%%% Cluster
				Z = linkage(Y, vLinkMethods{iLinkMethod});
				T = cluster(Z,'maxclust',nClusters);
				[T, nT] = sort_cluster_labels(T);		% Always sort so largest cluster first

				for i = 1:nGlyphs
					vGlyphPolygons{i,6} = T(i);
				end


	% 			% Show clustered glyphs in situ
	% 			figure(iFig);
	% 			clf; hold on;
	% 
	% 			vColors = lapham( nClusters );
	% 			for i = 1:nClusters
	% 					iCluster = i;%iSortedClusters'
	% 
	% 				iClusterGlyphs = find(T==iCluster)';
	% 
	% 				for j = iClusterGlyphs
	% 					fill(vGlyphPolygons{j,2}(:,1)+vGlyphPolygons{j,4}(1), -vGlyphPolygons{j,2}(:,2)-vGlyphPolygons{j,4}(2), vColors(i,:));
	% 					plot(vGlyphPolygons{j,3}(:,1)+vGlyphPolygons{j,4}(1), -vGlyphPolygons{j,3}(:,2)-vGlyphPolygons{j,4}(2), '-', 'Color', [1 1 1]/2);
	% 	% 					text(vGlyphPolygons{j,4}(2),-vGlyphPolygons{j,4}(1), 1000, sprintf('%d',i),'Color',vColors(i,:),'BackgroundColor',[1 1 1]*7/8)
	% 				end
	% 
	% 			end
	% 			axis equal;
	% 			zoom on;
	% 			hold off;
	% 			getframe;
	% 			% zoom(4);

	%%
				% Show clustered glyphs with first (two or three) dimensions
				if(any(iVectorType==[1 2])) % Only works when X is a set of vectors

					figure(iFig+100);
					clf;
					hold on;
					nDimPlot = 2;

					fMinDist = get_distance_metrics(X);
					iScale = fMinDist/iMaxDim*40;

					vColors = lapham( nClusters );
					for i = 1:nClusters
						iClusterGlyphs = find(T==i)';
						for j = iClusterGlyphs
							if(nDimPlot==2)
								fill(vGlyphPolygons{j,2}(:,1)*iScale+X(j,1), -vGlyphPolygons{j,2}(:,2)*iScale+X(j,2), vColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
								zoom on;
							elseif(nDimPlot==3)
								fill3d(vGlyphPolygons{j,3}(:,1)*iScale+X(j,1), -vGlyphPolygons{j,3}(:,2)*iScale+X(j,2),vGlyphPolygons{j,3}(:,1)*0+X(j,3), vColors(i,:), 'EdgeColor', 'none');
							end
						end
					end
					axis equal;
					title(sprintf('Demotic Classifier Clustered %id Vectors (Vector Type: %s, Linkage: %s)', nDimPlot, strVectorType, vLinkMethods{iLinkMethod}));
					hold off;
					getframe;

					savewipimage(sprintf('Demotic_Classifier_Clustered_Vectors (Vector Type - %s, Linkage - %i %s)', strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
				end
	% 
	% 			% Show glyphs overlaid in clusters
	% 			figure(iFig+200);
	% 			clf;
	% 
	% 			vColors = lapham( nClusters );
	% 			for i = 1:nClusters
	% 				subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
	% 				hold on;
	% 				iCluster = i;
	% 				iClusterGlyphs = find(T==iCluster)';
	% 				for j = iClusterGlyphs
	% 					plot(vGlyphPolygons{j,2}(1:end-1,1), -vGlyphPolygons{j,2}(1:end-1,2), '-', 'Color', vColors(i,:));
	% 				end
	% 				title(sprintf('%i',iCluster));
	% 				hold off;
	% 				axis off;
	% 				axis equal;
	% 			end
	% 			hold off;
	% 			sgtitle(sprintf('Demotic Classifier Clusters (Vector Type: %s, Linkage: %s)', strVectorType, vLinkMethods{iLinkMethod}));
	% 			getframe;
	% 
	% 			savewipimage(sprintf('Demotic_Classifier_Clusters (Vector Type - %s, Linkage - %i %s)', strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
			end
		end

	end

	%% Save working data for next step
	save('glyph_data','vGlyphs', 'vGlyphPolygons');

end

















