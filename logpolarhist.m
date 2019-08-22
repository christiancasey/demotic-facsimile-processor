function m = logpolarhist(vXY, nTh, nR, maxR)
% Returns a log-polar histogram of the shape based on a list of points in shape 
% nTh defines the number of angular bins, nR the number of magnitudes

	% Center the vectors first
	vC = mean(vXY,1);
	if(~all(vC==0))
		for i = 1:2
			vXY(:,i) = vXY(:,i) - vC(i);
		end
	end

	vTh = getangle(vXY);
	vR = log(getmagnitude(vXY)+1);
	
	if(~exist('nTh'))
		nTh = 16;
	end
	if(~exist('nR'))
		nR = 8;
	end	
	if(~exist('maxR'))
		maxR = max(vR);
	end
	
	vRBins = linspace(0, maxR, nR+1);
	vThBins = linspace(0,2*pi,nTh+1);
	
	m = zeros(nR,nTh);
	
	for i = 1:nR
		iR = (vRBins(i) < vR) & (vR <= vRBins(i+1));
		
		m(i,:) = histcounts(vTh(iR), vThBins);
		
	end













