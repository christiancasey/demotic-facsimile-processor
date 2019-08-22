function vCM = lapham(m, fSatScale)
%	LAPHAM
%   Based on IRIS, 
%	color scheme taken from visualizations in Lapham's Quarterly

	if nargin < 1
	   m = size(get(gcf,'colormap'),1);
	end
	if nargin < 2
	   fSatScale = 1.0;
	end 


	iTh = linspace(0,pi*2-0.1,m);
	vCM = (fSatScale*[ cos(iTh)' cos(iTh-2*pi/3)' cos(iTh-4*pi/3)' ] + 1 )/2;

	vCM = vCM * diag([0.5 0.5 0.5]) + ones(m,1)*[0.25 0.25 0.25];


end















