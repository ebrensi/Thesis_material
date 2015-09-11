function [lambda, W, rr_norm, resids] = band_Arnoldi_Ritz(result)
		% Exctract structures from this run of band-Arnoldi
		nn = length(result.H);
		N = size(result.V,1);
		H = result.H;
		Vh_defl = result.Vh_defl;
		Nn = size(Vh_defl,1);
		Iv = result.Iv;
		mc = result.mc;
		normH_est = result.tol.normA;
		
		% Construct candidate vector residual rVh = Vh*[0 0 ... 0 I]
		Vh = Vh_defl(:,Iv.ph); % extract remaining candidate vectors
		if nn > mc
			F2 = [zeros(mc,nn-mc) eye(mc)];
		end
		
		% construct deflated vector residual rVd = Vd*[0 ... ej ... ek 0...0],
		%   where j,k are locations of deflations.
		d = Iv.nd;  % number of deflated vectors
		Vd = zeros(N,0);
		F1 = zeros(d,nn);
		if d > 0
			Vd = Vh_defl(:,Iv.pd); % deflated vectors
			Ipos = Iv.I>0;  % positions of deflated vectors w/ positive index
			F1(Ipos, Iv.I(Ipos) ) = 1;
		end
		
		% Construct vector fT of total residual comlumn norms
		[W D] = eig(H);  % Eigenvalue decomposition of Rayleigh-quotient H
		lambda = diag(D);
		resids = (Vd*F1 + Vh*F2)*W;
		fT = colnorms(resids);
		rr_norm = fT ./ abs(lambda)';
	end
