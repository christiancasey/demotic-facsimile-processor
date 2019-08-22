

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
load glyph_data02;


%% Load parameter string lists
dimensionality_reduction_methods;		% Load list of DR methods: vDRMethods
distance_metrics_linkage_methods;		% Load vDistMetrics and vLinkMethods

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

			pause(10);
		end
	end

	save(sprintf('DimensionalityReduction_%03i',nDim),'vDRX');


end

%% Display Glyphs with Reduced Dimensions

manual_cluster_assignments; % Load variable with manual clusters
vManualClusterAssignments = vManualClusterAssignments(:,2);
vColors = lapham( max(vManualClusterAssignments) );

close all;
iFig = 1;

nDim = 2;
load(sprintf('DimensionalityReduction_%03i',nDim));

for iDRMethod = 1:34
	
	for iVectorType = 1:2
		
	
		if(iVectorType==1)
			X = mGlyphPolygons;
			strVectorType = 'Polygon';
		elseif(iVectorType==2)
			X = mGlyphHistograms;
			strVectorType = 'Histogram';
		end

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
		iScale = fMedianDist/get_glyph_extent(vGlyphs(:,5))*5;

		
		for j = 1:nGlyphs
			if(vManualClusterAssignments(j)>0)		% Don't draw clusters that haven't been assigned
				fill(vGlyphPolygons{j,2}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,2}(:,2)*iScale+Xdr(j,2), vColors(vManualClusterAssignments(j),:), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
			end
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


manual_cluster_assignments; % Load variable with manual clusters
vManualClusterAssignments = vManualClusterAssignments(:,2);
vManualClusterAssignments = vManualClusterAssignments+1;


nClusters = max(vManualClusterAssignments);

nDim = 2;
load(sprintf('DimensionalityReduction_%03i',nDim));

% Normalized mutual information for each setting (DR Types, Linkages, Vector Types)
vNMI = zeros(34, 7, 2);

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

% 		if ~(iVectorType== 3)
% 			try
% 				vNumClusters = 1:64;
% 				vClusterEval = evalclusters(Xdr,'linkage','Silhouette','KList', vNumClusters);
% 
% 				[fEvalMax,iEvalMax] = max(vClusterEval.CriterionValues);
% 				nBestCluster = vNumClusters(iEvalMax)
% 
% 				figure(21);
% 				clf; hold on;
% 				plot(vNumClusters, vClusterEval.CriterionValues);
% 				plot(nBestCluster, fEvalMax, '*r');
% 				title(sprintf('Demotic Cluster Evaluation (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
% 				hold off;
% 
% 				% print(imagefilename('Latin_Cluster_Eval'), '-depsc', '-painters');
% 				savewipimage(sprintf('Demotic_Classifier_Cluster_Eval (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
% 			catch
% 				warning(sprintf('Error in Demotic Cluster Evaluation (Dimensionality Reduction Method - %s, Vector Type - %s)', vDRMethods{iDRMethod}, strVectorType));
% 			end
% 		end

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
				

				%% Show clustered glyphs with first (two or three) dimensions
				if(any(iVectorType==[1 2])) % Only works when X is a set of vectors

					% Calculate normalized mutual information to compare clustering to manual clusters
					vNMI(iDRMethod,iLinkMethod,iVectorType) = normalized_mutual_information(vManualClusterAssignments, T)
% 					for i = 1:nGlyphs
% 						vGlyphPolygons{i,6} = T(i);
% 					end
% 
% 					figure(iFig+1000);
% 					clf;
% 					hold on;
% 					nDimPlot = 2;
% 
% 					fMinDist = get_distance_metrics(Xdr);
% 					iScale = fMinDist/iMaxDim*100;
% 
% 					vColors = lapham( nClusters );
% 					for i = 1:nClusters
% 						iClusterGlyphs = find(T==i)';
% 						for j = iClusterGlyphs
% 							if(nDimPlot==2)
% 								fill(vGlyphPolygons{j,2}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,2}(:,2)*iScale+Xdr(j,2), vColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% 								zoom on;
% 							elseif(nDimPlot==3)
% 								fill3d(vGlyphPolygons{j,3}(:,1)*iScale+Xdr(j,1), -vGlyphPolygons{j,3}(:,2)*iScale+Xdr(j,2),vGlyphPolygons{j,3}(:,1)*0+Xdr(j,3), vColors(i,:), 'EdgeColor', 'none');
% 							end
% 						end
% 					end
% 					axis equal;
% 					title(sprintf('Demotic Classifier Clustered %id Vectors (Dimensionality Reduction Method - %s, Vector Type: %s, Linkage: %s)', nDimPlot, vDRMethods{iDRMethod}, strVectorType, vLinkMethods{iLinkMethod}));
% 					hold off;
% 					getframe;
% 
% 					savewipimage(sprintf('Demotic_Classifier_Clustered_Vectors (Dimensionality Reduction Method - %s, Vector Type - %s, Linkage - %i %s)', vDRMethods{iDRMethod}, strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
% 				end
% 	%%
% 				% Show glyphs overlaid in clusters
% 				figure(iFig+200);
% 				clf;
% 	
% 				vColors = lapham( nClusters );
% 				for i = 1:nClusters
% 					subplot(ceil(sqrt(nClusters)),ceil(sqrt(nClusters)),i);
% 					hold on;
% 					iCluster = i;
% 					iClusterGlyphs = find(T==iCluster)';
% 					for j = iClusterGlyphs
% 						plot(vGlyphPolygons{j,2}(1:end-1,1), -vGlyphPolygons{j,2}(1:end-1,2), '-', 'Color', vColors(i,:));
% 					end
% 					title(sprintf('%i',iCluster));
% 					hold off;
% 					axis off;
% 					axis equal;
% 				end
% 				hold off;
% 				sgtitle(sprintf('Demotic Classifier Clusters (Dimensionality Reduction Method - %s, Vector Type: %s, Linkage: %s)', vDRMethods{iDRMethod}, strVectorType, vLinkMethods{iLinkMethod}));
% 				getframe;
% 	
% 				savewipimage(sprintf('Demotic_Classifier_Clusters (Dimensionality Reduction Method - %s, Vector Type - %s, Linkage - %i %s)', vDRMethods{iDRMethod}, strVectorType, iLinkMethod, vLinkMethods{iLinkMethod}));
			end
		end

	end

end



















