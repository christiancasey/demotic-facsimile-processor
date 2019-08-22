

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
% vManualClusterAssignments = vManualClusterAssignments(:,2);
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
		iScale = fMedianDist/get_glyph_extent(vGlyphs(:,6))*5;

		
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

% 		savewipimage(sprintf('Demotic Classifier %id Reduced Vectors (Dimensionality Reduction Method - %s -> %id vectors, Vector Type: %s)', nDimPlot, vDRMethods{iDRMethod}, nDim, strVectorType));	
	end
end
	
	
	
	
	
	
	
