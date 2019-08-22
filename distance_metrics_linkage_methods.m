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

vVectorTypes = { 'Polygon', 'Histogram' };