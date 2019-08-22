function [T, nT] = sort_cluster_labels(T)
% Rearranges cluster labels so that the largest cluster is number one, etc...

	nClusters = max(T);
	
	% Sort clusters according to number of members
	nT = histc(T,1:nClusters);
	[nT,iSortedClusters] = sort(nT,'descend');
	
	% Relabel clusters so that the most populous is number one, etc...
	T_sorted = T*0;
	for i = 1:nClusters
		T_sorted(T==iSortedClusters(i)) = i;
	end
	T = T_sorted;