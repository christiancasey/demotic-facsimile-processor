	function vCM = iris(m)
	%	IRIS
	%   Similar to Jet and HSV, but changes less abruptly, and so looks much less harsh
	%	Follows the order of the visual light spectrum like HSV, but more closely resembles a rainbow

		if nargin < 1
		   m = size(get(gcf,'colormap'),1);
		end

		fPurple = 2*pi/3;

		phi = 0;
		iTh = linspace(0,pi*2-fPurple,m);
		vCM = ([ cos(iTh+phi)' cos(iTh-2*pi/3+phi)' cos(iTh-4*pi/3+phi)' ] + 1 )/2;

	end















