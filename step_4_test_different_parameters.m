

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


%% Clustering
warning('off','all')

manual_cluster_assignments; % Load variable with manual clusters



nClusters = max(vManualClusterAssignments);


% Normalized mutual information for each setting
%	( Dimensions, Vector Type, Dimensionality Reduction Method, Distance Metric, Linkage )
vNMI = zeros(4,2,34,12,7);


vDims = [2 5 10 20];
for iDim = 1:4
	disp(iDim);
	
	nDim = vDims(iDim);
	load(sprintf('DimensionalityReduction_%03i',nDim));

	for iVectorType = 1:2

		strVectorType = vVectorTypes{iVectorType};

		for iDRMethod = 1:34

			Xdr = abs(vDRX{iDRMethod,iVectorType});

			% Skip DR methods that didn't work
			if(Xdr == 0)
				continue;
			end
			if( size(Xdr,2) < 2 )
				continue;
			end

			% Cluster
			for iDistMetric = 1:12
				disp(sprintf('%i %i %i %i', iDim, iVectorType,iDRMethod,iDistMetric));

				try
					Y = pdist(Xdr,vDistMetrics{iDistMetric});
				catch
					continue;
				end

				for iLinkMethod = 1:7

% 					disp(sprintf('%i %i %i %i %i', iDim, iVectorType,iDRMethod,iDistMetric,iLinkMethod));
					% Cluster
					Z = linkage(Y, vLinkMethods{iLinkMethod});
					T = cluster(Z,'maxclust',nClusters);
					[T, nT] = sort_cluster_labels(T);		% Always sort so largest cluster first

					% Calculate normalized mutual information to compare clustering to manual clusters
					vNMI(iDim, iVectorType,iDRMethod,iDistMetric,iLinkMethod) = normalized_mutual_information(vManualClusterAssignments, T);

					
				end

			end
		end

	end

end

%% Output best n results
[vMaxNMI, iMaxNMI] = sort(vNMI(:),'descend');

vBest = zeros(length(iMaxNMI),5);
vBestLabels = cell(length(iMaxNMI),7);

for i = 1:length(iMaxNMI)
	[iDim,iVT,iDR,iDM,iLM] = ind2sub(size(vNMI),iMaxNMI(i));
	vBest(i,:) = [iDim,iVT,iDR,iDM,iLM];
	vBestLabels(i,:) = { i, vMaxNMI(i), vDims(vBest(i,1)), vVectorTypes{vBest(i,2)}, vDRMethods{vBest(i,3)}, vDistMetrics{vBest(i,4)}, vLinkMethods{vBest(i,5)} };
end

tBestLabels = cell2table(vBestLabels);%,'VariableNames',{'Rank', 'NMI', 'No. of Dimensions', 'Vector Type', 'Dim. Red. Method', 'Distance Metric', 'Linkage'});
writetable(tBestLabels,'Best Clusterings.csv');

%% Show clusters for best methods

n = 1;
for i = 1:n

	nDim = vDims(vBest(i,1));
	load(sprintf('DimensionalityReduction_%03i',nDim));
	strVectorType = vVectorTypes{vBest(i,2)};
	Xdr = abs(vDRX{vBest(i,3),vBest(i,2)});

	% Skip DR methods that didn't work
	if(Xdr == 0)
		error('Invalid DR Method');
	end
	if( size(Xdr,2) < 2 )
		error('Invalid DR Method');
	end

	try
		Y = pdist(Xdr,vDistMetrics{vBest(i,4)});
	catch
		error('Invalid Distance Metric');
	end


	% Cluster
	Z = linkage(Y, vLinkMethods{vBest(i,5)});
	T = cluster(Z,'maxclust',nClusters);
	[T, nT] = sort_cluster_labels(T);		% Always sort so largest cluster first

	[ vMaxNMI(i) normalized_mutual_information(vManualClusterAssignments, T) ]

	[fMinDist, fMeanDist, fMedianDist, fMaxDist ] = get_distance_metrics(Xdr);
	iScale = fMedianDist/get_glyph_extent(vGlyphs(:,6))*1;

	figure(i);

% 	plot_glyph_scatter(Xdr, vGlyphPolygons(:,2), T);
% 	sgtitle(sprintf('Demotic Classifier Best Parameters, Rank %i', i));

	plot_clusters(vGlyphPolygons(:,2), T);
	sgtitle(sprintf('Demotic Classifier Best Parameters, Rank %i', i));
	
	zoom on;
	getframe;

% 	savewipimage(sprintf('Demotic Classifier Best Parameters, Rank %i', i));	
	

end


%% Create view of manual assignments for comparison

figure(101);
vManualClusterAssignments = sort_cluster_labels(vManualClusterAssignments);
plot_clusters(vGlyphPolygons(:,2), vManualClusterAssignments);
sgtitle('Demotic Classifier Manual Cluster Assignments');
% savewipimage('Demotic Classifier Manual Cluster Assignments');


%% Create sign list



vGlyphLocations = cell2mat(vGlyphs(:,1:5));




n = 1;
for i = 1:n

	nDim = vDims(vBest(i,1));
	load(sprintf('DimensionalityReduction_%03i',nDim));
	strVectorType = vVectorTypes{vBest(i,2)};
	Xdr = abs(vDRX{vBest(i,3),vBest(i,2)});

	% Skip DR methods that didn't work
	if(Xdr == 0)
		error('Invalid DR Method');
	end
	if( size(Xdr,2) < 2 )
		error('Invalid DR Method');
	end

	try
		Y = pdist(Xdr,vDistMetrics{vBest(i,4)});
	catch
		error('Invalid Distance Metric');
	end


	% Cluster
	Z = linkage(Y, vLinkMethods{vBest(i,5)});
	T = cluster(Z,'maxclust',nClusters);
	[T, nT] = sort_cluster_labels(T);		% Always sort so largest cluster first
	
	strFilename = sprintf('Classifier Sign List, Rank %i.tex',i);
	output_sign_list(strFilename, vGlyphLocations, vGlyphPolygons(:,2), T);
end




































































































	
	