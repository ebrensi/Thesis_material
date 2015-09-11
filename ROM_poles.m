	function [mu wt] = ROM_poles(cV,W,rho,lambda,s0,FRdomain)
		mu = 1./lambda + s0;
		% Compute pole mass
		f = colnorms(cV*W)';
		g = colnorms((W\rho).')';
		delta = abs(s0-mu) ./ (d2S(mu,FRdomain)+1);% note this
		delta(isinf(mu)) = 1;
		wt = delta .* f .* g;
	end
