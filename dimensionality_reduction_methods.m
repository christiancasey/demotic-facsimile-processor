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