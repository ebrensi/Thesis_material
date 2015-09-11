function  [V ROM_realization] = ROM_bArnoldi(URM_realization, S0, J)

% parameters
dtol = sqrt(eps);
ctol = dtol;
keep_tol = ctol;
rel_wt_tol = inf;
ROI = [9 10];

% Model setup
A = ROM_realization.A;
E = ROM_realization.E;
B = ROM_realization.B;
C = ROM_realization.C;

N = size(A,1);
VROM =  zeros(N,0);
Y = zeros(N,0);
flops = 0;

for j = 1:length(S0)
	s0 = S0(j);
	ell(j) = size(Y,2);
	[multH R] = make_SI_op(A,E,B,s0);
	m = size(R,2) + ell(j);
	fprintf('\ncycle %d expanding at s_0 = %s,  band_size = %d + %d \n',j,s0string(s0),size(R,2),ell(j));
	
	if isempty(Y)
		result = band_Arnoldi(m,R,multH);
	else
		% Process Y (thick-starting basis) and manually deflate the resulting candidates
		result = band_Arnoldi(ell(j),[Y R],multH);
		result.Vh_defl(:,1:ell(j)) = 0;  % force deflation of thick-started Ritz-vectors
	end
	
	% Continue iterating to J(j)
	n0 = size(result.V,2);
	result = band_Arnoldi(J(j)+ell(j),[Y R],multH,[],n0,result);
	n = size(result.V,2);
	
	% count flops
	if isreal(s0)
		flops = flops + result.flops;
	else
		flops = flops + 4*result.flops;
	end
	
	
	% Determine Ritz-values for this cycle
	[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
	rho2 = result.rho(:,(ell(j)+1):m);
	cV = C'*result.V;
	[mu wt] = ROM_poles(cV, W, rho2, lambda, s0, FRdomain);
	
	
	converged = rr' < ctol;
	ninf = d2S(mu,s) < 1e10;
	in_ROI = imag(mu) > 0 & imag(mu) < 1e11;
	
	% Set criteria for Ritz vectors to keep
	rel_wt = wt / sum(wt);
	keep = (rr' < keep_tol | rel_wt > rel_wt_tol);
	
	% 	[rrsort idx] = sort(rr);
 	% 	fprintf('&$\\Re(\\mu)$ & $\\Im(\\mu)$ &\\texttt{rr} & wt & keep \\\\ \n \\hline \n')
	% 	for i = 1:min(10,length(mu))
	% 		fprintf('%d & \\texttt{%4.4e} &  \\texttt{%4.4e}i & %g & %g & %d \\\\\n',...
	% 			i, real(mu(idx(i))),imag(mu(idx(i))),rr(idx(i)), wt(idx(i)), keep(idx(i)) );
	% 	end
	% 	fprintf('\\hline\n');
	ntot = size(VROM,2)+n-ell(j);
	fprintf('... ROM: %d, n_tot: %d,   converged: %d,   keep: %d  weight: %g\n',...
		n,ntot,nnz(converged), nnz(keep), sum(wt) );
	
	fprintf('...updating thick-restart basis...')
	Znew = result.V*W(:,keep);
	% 	Znew = [result.V*W(:,~converged & keep) real(result.V*W(:,converged & keep)) imag(result.V*W(:,converged & keep))];
	[Y T] = add2basis(Y, T, Znew,dtol);
	fprintf('dim Y = %d\n',size(Y,2));
	
	
	VROM = [VROM result.V(:,ell(j)+1:n)];
	j = j+1;
end % of j-th cycle of outer loop


%% Post basis-building analysis
[V VROM_split] = make_basis_real(VROM);
nreal = size(V,2);


%  Explicitly project system realization (A,E,B,C) onto V
ROM_realization.A = V'*A*V;
ROM_realization.E = V'*E*V;
ROM_realization.B = V'*B;
ROM_realization.C = V'*C;

end % main

