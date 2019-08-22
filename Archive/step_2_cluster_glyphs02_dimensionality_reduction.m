

%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap lapham;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');

iFig = 0;

%% Compare manual counts

vCounts = [ 45 162 73 17 12 81 3 6 1 42 18 28 13 9 8 6 27 3 21 2 3 3 10 1 10 11 2 1 1 1 4 4 3 1 5 3 2 6 1 1 1 1 1 1 1 1 1 3 1 1 1 10 1 1 2 1 3 1 ];
vCounts = vCounts';
vCounts = sort(vCounts, 'descend');
sum(vCounts)

%% Load parameter string lists
dimensionality_reduction_methods;		% Load list of DR methods: vDRMethods
distance_metrics_linkage_methods;		% Load vDistMetrics and vLinkMethods

%% Load working data from previous step

load glyph_data;

nGlyphs = size(vGlyphs,1);


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


%% Try clustering with different types of vectors to see which one works best


for nDim = [2 5 10 20]

	vDRX = cell(34,2);
	for iVectorType = 1:2

		if(iVectorType==1)
			X = mGlyphPolygons;
			strVectorType = 'Polygon';
		elseif(iVectorType==2)
			X = mGlyphHistograms;
			strVectorType = 'Histogram';
		end

		%% Dimensionality Reduction	
		for iDRMethod = 1:34

			disp(vDRMethods{iDRMethod});

			try
				[Xdr, mapping] = compute_mapping(X, vDRMethods{iDRMethod}, nDim);
				vDRX{iDRMethod,iVectorType} = Xdr;
			catch
				warning(sprintf('Error with method %i - %s', iDRMethod, vDRMethods{iDRMethod}));
				vDRX{iDRMethod,iVectorType} = 0;
				continue;
			end

			pause(60);
		end
	end

	save(sprintf('DimensionalityReduction_%03i',nDim),'vDRX');


end

return;

%% Display Glyphs with Reduced Dimensions
close all;
iFig = 1;

nDim = 2;
load(sprintf('DimensionalityReduction_%03i',nDim));

for iVectorType = 1:2
	
	if(iVectorType==1)
		X = mGlyphPolygons;
		strVectorType = 'Polygon';
	elseif(iVectorType==2)
		X = mGlyphHistograms;
		strVectorType = 'Histogram';
	end
	
	for iDRMethod = 1:34

		Xdr = abs(vDRX{iDRMethod,iVectorType});
		
		% Skip DR methods that didn't work
		if(Xdr == 0)
			continue;
		end
		if( size(Xdr,2) < 2 )
			continue;
		end
		
		%% Show clustered glyphs with first (two or three) dimensions
		iFig = iFig+1;
		figure(iFig+100);
		clf;
		hold on;
		nDimPlot = 2;

		[fMinDist, fMeanDist, fMedianDist, fMaxDist ] = get_distance_metrics(Xdr);
		iScale = fMedianDist/iMaxDim*5;

		vColors = lapham( nGlyphs );
		for j = 1:nGlyphs
			fill(vGlyphPolygons{j,2}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,2}(:,2)*iScale+Xdr(j,2), vColors(j,:), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
		end
		axis equal;
		axis off;
		title(sprintf('Demotic Classifier %id Reduced Vectors (Dimensionality Reduction Method - %s -> %id vectors, Vector Type: %s)', nDimPlot, vDRMethods{iDRMethod}, nDim, strVectorType));
		hold off;
		zoom on;
		getframe;

		savewipimage(sprintf('Demotic Classifier %id Reduced Vectors (Dimensionality Reduction Method - %s -> %id vectors, Vector Type: %s)', nDimPlot, vDRMethods{iDRMethod}, nDim, strVectorType));
				
	end
end
		

%% Clustering

for iVectorType = 1:2
	
	if(iVectorType==1)
		X = mGlyphPolygons;
		strVectorType = 'Polygon';
	elseif(iVectorType==2)
		X = mGlyphHistograms;
		strVectorType = 'Histogram';
	end
	
	for iDRMethod = 1:34

		Xdr = abs(vDRX{iDRMethod,iVectorType});
		
		% Skip DR methods that didn't work
		if(Xdr == 0)
			continue;
		end
		if( size(Xdr,2) < 2 )
			continue;
		end
		
	
		%% Evaluate clustering algorithms
		nClusters = 64;

		if ~(iVectorType== 3)
			try
				vNumClusters = 1:64;
				vClusterEval = evalclusters(Xdr,'linkage','Silhouette','KList', vNumClusters);

				[fEvalMax,iEvalMax] = max(vClusterEval.CriterionValues);
				nBestCluster = vNumClusters(iEvalMax)

				figure(21);
				clf; hold on;
				plot(vNumClusters, vClusterEval.CriterionValues);
				plot(nBestCluster, fEvalMax, '*r');
				title(sprintf('Demotic Cluster Evaluation (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
				hold off;

				% print(imagefilename('Latin_Cluster_Eval'), '-depsc', '-painters');
				savewipimage(sprintf('Demotic_Classifier_Cluster_Eval (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
			catch
				warning(sprintf('Error in Demotic Cluster Evaluation (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
			end
		end

		%% Cluster
		for iDistMetric = 2

			if(any(iVectorType==[1 2]))
				Y = pdist(Xdr,vDistMetrics{iDistMetric});
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

				%% Show clustered glyphs with first (two or three) dimensions
				if(any(iVectorType==[1 2])) % Only works when X is a set of vectors

					figure(iFig+1000);
					clf;
					hold on;
					nDimPlot = 2;

					fMinDist = get_distance_metrics(Xdr);
					iScale = fMinDist/iMaxDim*100;

					vColors = lapham( nClusters );
					for i = 1:nClusters
						iClusterGlyphs = find(T==i)';
						for j = iClusterGlyphs
							if(nDimPlot==2)
								fill(vGlyphPolygons{j,2}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,2}(:,2)*iScale+Xdr(j,2), vColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
								zoom on;
							elseif(nDimPlot==3)
								fill3d(vGlyphPolygons{j,3}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,3}(:,2)*iScale+Xdr(j,2),vGlyphPolygons{j,3}(:,1)*0+Xdr(j,3), vColors(i,:), 'EdgeColor', 'none');
							end
						end
					end
					axis equal;
					title(sprintf('Demotic Classifier Clustered %id Vectors (Dimensionality Reduction Method - %s, Vector Type: %s, Linkage: %s)', nDimPlot, vDRMethods{iDRMethod}, strVectorType, vLinkMethods{iLinkMethod}));
					hold off;
					getframe;

					savewipimage(sprintf('Demotic_Classifier_Clustered_Vectors (Dimensionality Reduction Method - %s, Vector Type - %s, Linkage - %i %s)', vDRMethods{iDRMethod}, strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
				end
	%%
				% Show glyphs overlaid in clusters
				figure(iFig+200);
				clf;
	
				vColors = lapham( nClusters );
				for i = 1:nClusters
					subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
					hold on;
					iCluster = i;
					iClusterGlyphs = find(T==iCluster)';
					for j = iClusterGlyphs
						plot(vGlyphPolygons{j,2}(1:end-1,1), -vGlyphPolygons{j,2}(1:end-1,2), '-', 'Color', vColors(i,:));
					end
					title(sprintf('%i',iCluster));
					hold off;
					axis off;
					axis equal;
				end
				hold off;
				sgtitle(sprintf('Demotic Classifier Clusters (Dimensionality Reduction Method - %s, Vector Type: %s, Linkage: %s)', vDRMethods{iDRMethod}, strVectorType, vLinkMethods{iLinkMethod}));
				getframe;
	
				savewipimage(sprintf('Demotic_Classifier_Clusters (Dimensionality Reduction Method - %s, Vector Type - %s, Linkage - %i %s)', vDRMethods{iDRMethod}, strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
			end
		end

	end

end

















