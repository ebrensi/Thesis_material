function  test

% parameters
data_filename = 'ex308'; 
S0 =  5i*2*pi*1e9;
dtol = sqrt(eps);
ctol = 1e-9;
tfunc_tol = 1e-2;
ROI = [9 10];
nmax = 100;
jmax = 10;

% Model setup
model_name = sprintf('%s',data_filename);
l = length(S0);
[s FRdomain] = tfdomain(200,ROI);
URM_FR = URM_freq_response(data_filename,[],[],FRdomain);
[A E C B] = realization(data_filename);
N = size(A,1);
cV = zeros(size(C,2),nmax);
Vreal = zeros(N,nmax*2);

VROM =  zeros(N,0);
kept_Ritz_vals = zeros(0,1);
kept_Ritz_vecs = zeros(N,0);
j = 1;
while j <= length(S0) && j < jmax  
	s0 = S0(j);
	ell = size(kept_Ritz_vecs,2);
	[multH R] = make_SI_op(A,E,B,s0);
	m = size(R,2) + ell;
	fprintf('%s: Band-Arnoldi cycle %d expanding at s_0 = %s,  m=%d\n',model_name,j,s0string(s0),m);
	
	result = band_Arnoldi(m,[kept_Ritz_vecs R],multH);  % Process the initial band
	n = size(result.V,2);
% 	result.Vh_defl(:,1:ell) = 0;
	cV(:,1:n) = C'*result.V;
	converged = [];
	ROM_err = inf;
	
	while nnz(converged) < ell+1  && n < nmax % inner loop
		% perform the (n+1)-th iteration 
		result = band_Arnoldi(n+1,[kept_Ritz_vecs R],multH,[],n,result);
		if length(result.H) ~= n+1
			break;
		end
		n = n+1;
		vn = result.V(:,n); 
		cV(:,n) = C'*vn;

		
		%% **** Inner-loop anaylsis *****
		%% Compute implicit ROM tfunc poles/weights
		[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
		converged = rr < ctol;
		[mu wt] = ROM_poles(cV(:,1:n), W, result.rho, lambda, s0, FRdomain);
		ROM_wt(n) = sum(wt);
	end
	
	% This code executes after every cycle at one point s0(j)
	[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
	[mu wt] = ROM_poles(cV(:,1:n), W, result.rho, lambda, s0, FRdomain);
	ROM_wt(n) = sum(wt);
	converged = rr' < ctol;
	dom = wt/sum(wt) > 1e-2;
	inf_= abs(lambda)/norm(result.H) < 1e-3;
	plottable = ~inf_; 
	keep = converged;
	normH_est = result.tol.normA;
	
	%% determine next interpolation point
	sig = ~converged & ~inf_; 
	mu_sig = mu(sig);
	[max_wt idx] = max(wt(sig));
	S0(j+1) = 1i*imag(mu_sig(idx));
 	
	
	%% Plot finite ROM poles
	plot_ROMpoles(mu(plottable), s0, rr(plottable), wt(plottable), ROI);
	hold on
	plot(real(S0(j+1)),imag(S0(j+1)),'ro');
	hold off
	title(sprintf('dominant poles,s_0=%s, n=%d, j=%d',s0string(s0),n,j));
	colorbar;	
	

	Z = result.V*W;
	kept_Ritz_vecs = Z(:,keep);
% 	kept_Ritz_vecs = [real(Z(:,keep)) imag(Z(:,keep))];

	
	fprintf('\t n=%d:  %d Ritz-pairs to keep:\n',n,nnz(keep));
	disp(sortrows([mu(keep) rr(keep)' wt(keep)]));
	
%  	H1_residuals = multH(Z) - Z*diag(lambda);
% 	AE_residuals = A*Z - E*Z*diag(mu);
% 	H1_residuals_bA = Vh*F2*W;
% 	
%	H1_rnorm = colnorms(H1_residuals)' ./ abs(lambda);
% 	H1_rnorm_bA = colnorms(H1_residuals_bA)' ./ abs(lambda);
% 	H1_rnorm_bA_est1 = norm(Vh,'fro')*colnorms(F2*W)' ./ abs(lambda);
% 	H1_rnorm_bA_est2 = (colnorms(Vh)*F2*abs(W))' ./ abs(lambda);
% 	AE_rnorm = colnorms(AE_residuals)' ./ colnorms(A*Z)';
% 	AE_rnorm_bA_est1 = colnorms(H1_residuals_bA)' ./(abs(mu).*abs(lambda)*norm(result.H));
% 	AE_rnorm_bA_est2 = norm(A-s0*E,'fro')*colnorms(H1_residuals_bA)' ./ (norm(A,'fro')*abs(lambda));
% 	figure;
% 	semilogy([H1_rnorm H1_rnorm_bA...
% 		H1_rnorm_bA_est1, H1_rnorm_bA_est2,...
% 		AE_rnorm AE_rnorm_bA_est1 AE_rnorm_bA_est2],'-o');
% 	line([1 length(lambda)],[dtol dtol],'LineStyle',':');
% 	title(sprintf('%s, s0=%s, n=%d:  rel-residual norms of some Ritz vectors',model_name,s0string(s0),n));
% 	legend('H_1 exp','H_1 bArnoldi',...
% 		'H_1 bA est1','H_1 bA est2',...
% 		'(A,E) exp','(A,E) bA est1','(A,E) bA est2');
% 	plottools('on','plotbrowser')
% 		
		
% 	if j < length(S0)
% 		s1 = S0(j+1);
% 		[multH2 R2] = make_SI_op(A,E,B,s1);
% 		Delta = (s1-s0);
% 		T = @(v) Delta*multH2(v) + v;
% 		zeta = 1./ ( 1 - Delta*lambda );
% 		lambda2 = zeta.*lambda ;
% 		
% 		% Compute actual residuals of kept Ritz-vectors
% % 		TH1_residuals = T(H1_residuals)*diag(zeta);
% % 		TH1_residuals_bA = T(H1_residuals_bA)*diag(zeta);
% % 		H2_residuals = multH2(Z) - Z*diag(lambda2);
% % 		H2_rnorm = colnorms(H2_residuals)' ./ abs(lambda2);
% % 		H2_rnorm_bA = colnorms(TH1_residuals_bA)' ./ abs(lambda2);
% % 		
% % 		figure;
% % 		semilogy([H1_rnorm H1_rnorm_bA...
% % 			H2_rnorm H2_rnorm_bA],'-o');
% % 		line([1 length(lambda)],[dtol dtol],'LineStyle',':');
% % 		title(sprintf('%s, s0=%s, n=%d:  rel-residual norms of some Ritz vectors',model_name,s0string(s0),n));
% % 		legend('H_1 exp','H_1 bArnoldi',...
% % 			'H_2 exp','H_2 bArnoldi');
% % 		plottools('on','plotbrowser')
% 	end

	
	%%  ******** Analysis of Intermediate ROM after j-th cycle  ***********
	
	% Compute explicit ROM frequency response and relative-error
	realV = make_basis_real(result.V);
	ROM_exp_FR = transfer_function(realV,A,E,C,B,s);
	tf_err = tfunc_err(URM_FR,ROM_exp_FR);
	plot_MIMO_tfunc(ROM_exp_FR, URM_FR,FRdomain,[1,1]);
    suplabel(sprintf('%s, n=%d: FR, err=%g,  s0 = %s',model_name,size(realV,2),tf_err,s0string(s0)),'t');

	VROM = [VROM result.V(:,ell+1:n)];
	j = j+1;
end % of j-th cycle of outer loop 

%% Determine ROM error for intermediate explicitly projected model
Vreal = [];
for i = 1:n		
	vn = result.V(:,n);
	Vreal = orth([Vreal real(vn) imag(vn)]);
	ROM_exp_FR = transfer_function(Vreal,A,E,C,B,s);
	[tfexp_err(n) fullerr] = tfunc_err(URM_FR,ROM_exp_FR);
	ROM_err = tfexp_err(n);
end
% Post basis-building analysis
figure
h = subplot(2,1,1);
semilogy(1:n,tfexp_err(1:n),'k.')
line([1 n],[tfunc_tol tfunc_tol],'LineStyle',':');		
title('rel-error');
set(h,'XTick',[]);

subplot(2,1,2);
semilogy(ROM_wt,'b.');
title('ROM weight');



	function [mu wt] = ROM_poles(cV,W,rho,lambda,s0,FRdomain)
		mu = 1./lambda + s0;
		% Compute pole mass
		f =  sum(abs(cV*W),1)';
		g = sum(abs(W\rho),2);
		delta = abs(s0-mu) ./ (d2S(mu,FRdomain)+1);% note this
		delta(isinf(mu)) = 1;
		wt = delta .* f .* g;
	end

	function [lambda, W, rr_norm, Vd,F1, Vh, F2] = band_Arnoldi_Ritz(result)
		% Exctract structures from this run of band-Arnoldi
		nn = length(result.H);
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
		Vd = [];
		F1 = zeros(d,nn);
		if d > 0
			Vd = Vh_defl(:,Iv.pd); % deflated vectors
			Ipos = Iv.I>0;  % positions of deflated vectors w/ positive index
			F1(Ipos, Iv.I(Ipos) ) = 1;
		end
		
		% Construct vector fT of total residual comlumn norms
		fT = colnorms(Vd*F1) + colnorms(Vh*F2);
		
		[W D] = eig(H);  % Eigenvalue decomposition of Rayleigh-quotient H
		lambda = diag(D);
		rr_norm = fT*abs(W) ./ abs(lambda)';
	end

	function [err fullerr] = tfunc_err(URM_FR,ROM_FR)
		Ls = size(URM_FR,3);
		fullerr = zeros(1,Ls);
		for i = 1:Ls
			fullerr(i) = norm(URM_FR(:,:,i) - ROM_FR(:,:,i)) / norm(URM_FR(:,:,i));
		end
		err = norm(fullerr);
	end
	function plot_MIMO_tfunc(ROM_tfunc, URM_tfunc,FRdomain,ROMs)
		nin = size(ROM_tfunc,1);
		Ls = size(ROM_tfunc,3);
		
		if ~exist('ROMs','var')
			plot_number = 1;
			tot = nin*(nin+1)/2;
			for i = 1:nin
				for k = 1:i
					URM = squeeze(URM_tfunc(i,k,:));
					ROM = squeeze(ROM_tfunc(i,k,:));
					h = subplot(tot,1,plot_number);
					semilogy(FRdomain,abs(URM),'r-+',FRdomain,abs(ROM));
					title(sprintf('(%d,%d), err=%g',i,k,norm(URM-ROM)/norm(URM)));
% 					saveas(gcf,sprintf('%s_%d%d.png',model_name,i,k))
					if plot_number < tot
						set(h,'XTick',[]);
					end
					plot_number = plot_number + 1;
				end
			end
		else
			tot = size(ROMs,1);
			for plot_number = 1:tot
				figure;
				i = ROMs(plot_number,1);
				k = ROMs(plot_number,2);
				URM = squeeze(abs(URM_tfunc(i,k,:)));
				ROM = squeeze(abs(ROM_tfunc(i,k,:)));
				semilogy(FRdomain,URM,'r--',FRdomain,ROM);
				title(sprintf('(%d,%d), err=%g',i,k,norm(URM-ROM)/norm(URM)));
			end
		end
	end

end % main

