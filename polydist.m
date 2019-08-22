function mDist = polydist(v)

	n = length(v);
	
	mDist = zeros(1,(n*(n-1))/2);
	
	k = 0;
	for i = 1:n
		for j = i+1:n
			k = k+1;
			mDist(k) = polydiff(v{i},v{j});
		end
	end