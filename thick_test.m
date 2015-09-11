function [H_full H_fd H_nd]  = thick_test

N = 100;
m = 2;
k = 15;
n_keep = 3; 

force_deflation = false;

dtol = sqrt(eps);

A = rand(N);
R = rand(N,m);
multA = @(v) A*v;

Y = [];  T = [];  kept_lambda = [];

for j = 1:3
	ell = size(Y,2);
	fprintf('N = %d, bandsize: %d + %d,   k=%d\n', N,ell, m, k);
	if j == 1
		result = band_Arnoldi(k,R,multA);
	else		
		result = band_Arnoldi(ell,[Y R],multA);		
		if force_deflation
			result.Vh_defl(:,1:ell) = 0;  % force deflation of thick-started Ritz-vectors
		end
		n0 = size(result.V,2);
		result = band_Arnoldi(k,[Y R],multA,[],n0,result);
	end
	H = result.H;
	abs(H)
	spyH(H)
	
	% Get Ritz-pairs
	[lambda, W, rr_norm, Vd,F1, Vh, F2] = band_Arnoldi_Ritz(result);
	[rr_sorted idx] = sort(rr_norm);
	keep = idx(1:n_keep);
	
	Z2keep = result.V * W(:,keep);
	lambda2keep = lambda(keep);
	nY_old = size(Y,2); 
	
	[Y T kept] = add2basis(Y, T, Z2keep,1e-4);
	kept_lambda = [kept_lambda lambda2keep(kept).'];
	nY_new = size(Y,2); 
	
	U = T * diag(kept_lambda)/T;


	fprintf('\nAdding %d Ritz-pairs for total of %d:\n',nY_new-nY_old, nY_new);
	disp([abs(kept_lambda)]);
end


	
end
