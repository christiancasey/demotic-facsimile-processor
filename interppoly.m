function vPi = interppoly(vP, nPts, cMethod)		
% Returns a polygon interpolated to a set number of points

	if(~exist('cMethod'))
		cMethod = 's';
	end
	
	vPx = vP(:,1);
	vPy = vP(:,2);
	[vPx, vPy] = polysplit(vP(:,1), vP(:,2));
	nP = size(vPx,1);
	if( nP <= 0 )
		vPi = [];
		return;
	elseif( nP == 1 )
		vP = [vPx{1}, vPy{1}];
	else
		sMag = sum(mag(vP(any(~isnan(vP),2),:) ));
		vPi = [];
		for j = 1:nP
			vPij = [vPx{j}, vPy{j}];
			nPi = round( sMag * nPts / sum(mag(vPij)) );
			nPi = max(nPi,2);
			nPi = min(nPi,nPts);
		
			vPij = interppoly(vPij, nPi, cMethod);
			if(j>1)
				vPi = [ vPi ; NaN NaN ];
			end
			vPi = [ vPi ; vPij ];
		end
		return;
	end
		
	[n d] = size(vP);
	ip = 1:n;
	ii = linspace(1,n,nPts);

	vPi = zeros(nPts,d);
	
	for i = 1:d

		vPd = vP(:,i);

		vPdI = interp1(ip,cumsum(vPd),ii,cMethod);
		
		vPdI = diff(vPdI)*nPts/n;
		vPdI = [ vPdI(end) ; vPdI(:) ];
		vPdI = vPdI - (mean(vPdI)-mean(vPd));
		
		vPi(:,i) = vPdI;
	end



return;
% 
% function vPi = interppoly(vP, nPts, cMethod)		
% 	
% 	if(~exist('cMethod'))
% 		cMethod = 's';
% 	end
% 	
% 	if( all(isnan(vP(end,:)),2) )
% 		vP = vP(1:end-1,:);
% 	end
% 	
% 	vPx = vP(:,1);
% 	vPy = vP(:,2);
% 	if( any( all(isnan(vP),2) ) )
% 		[vPx, vPy] = polysplit(vPx, vPy);
% 	else
% 		vPx = {vPx};
% 		vPy = {vPy};
% 	end
% 	
% 	sLen = sum(mag(vP));
% 	vP = [];
% 	
% 	nShapes = size(vPx,1);
% 	for i = 1:nShapes
% 		vPi = [vPx{i} vPy{i}];
% 		
% 		nPi = round( sLen * nPts / sum(mag(vPi)) );
% 		nPi = max(nPi,2);
% 		nPi = min(nPi,nPts);
% 		
% 		[n d] = size(vPi);
% 		ip = 1:n;
% 		ii = linspace(1,n,nPi);
% 
% 		vPij = zeros(nPi,d);
% 
% 		for j = 1:d
% 
% % 			vPd = vPi(:,j);
% 
% 			vPdJ = interp1(ip,cumsum(vPi(:,j)),ii,cMethod);
% 
% 			vPdJ = diff(vPdJ)*nPts/n;
% 			vPdJ = [ vPdJ(end) ; vPdJ(:) ];
% 			vPdJ = vPdJ - (mean(vPdJ)-mean(vPi(:,j)));
% 
% 			vPij(:,j) = vPdJ;
% 		end
% 		
% 		
% 		
% 		if(i>1)
% 			vP = [ vP ; NaN NaN ];
% 		end
% 		vP = [ vP ; vPij ];
% 	end
% 
% 
% 
% 
% 
% 
% 
