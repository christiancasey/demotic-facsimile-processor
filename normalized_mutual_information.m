function vNMI = normalized_mutual_information(c, v)
% Returns normalized mutual information between manually-assigned classes: c( No. of Objects, 1 )
%	and each in a set of automatically generated cluster labels: v( No. of Objects, No. of Clusterings )

% c = vManualClusterAssignments;
% v = T;

	n = length(c);

	if(n ~= size(v,1))
		error('Classes and clusters must comprise the same number of objects.');
	end
	
	m = size(v,2);						% No. of clusterings to test
	vNMI = zeros(1,m);					% Output variable
	
	nC = max(c);						% No. of classes
	pC = histc(c,1:nC) / n;				% Probability of classes
	hC = ( -1 * pC .* log(pC) );		% Entropy of classes
	hC(pC==0) = 0;							% Deal with log(0) problem
	hC = sum(hC);
	
	for i = 1:m
		
		nV = max(v(:,i));					% No. of clusters
		pV = histc(v(:,i),1:nV) / n;		% Probability of clusters
		hV = ( -1 * pV .* log(pV) );		% Entropy of clusters
		hV(pV==0) = 0;							% Deal with log(0) problem
		hV = sum(hV);
		
		pVC = zeros(nV,nC);					% Matrix for joint probilities
		for j = 1:nC
			for k = 1:nV
				pVC(k,j) = sum( (c==j) & (v(:,i)==k) ) / n;
			end
		end
		
		iVC = pVC .* log( pVC ./ (pV*pC') );	% Joint information
		iVC(pVC==0) = 0;						% Deal with log(0) problem
		
		vNMI(i) = 2 * sum(iVC(:)) / ( hC + hV );
	end
	
	% If there is only one clustering to test, output a single value
	if(m==1)
		vNMI = vNMI(1);
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	