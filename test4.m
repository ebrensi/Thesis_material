function  [ell flops V] = test4(data_filename, S0, J)
% [ell flops V] = test4(data_filename, S0, J)
%
% test4 takes preset interpolation point and # of iterations
%   it also performs comparisons with URM

% Parameters
dtol = sqrt(eps);
ctol = dtol;
keep_tol = 1e-4;
rel_wt_tol = inf;
ROI = [9 10]; axis3d = [0 10.1 9 10];  % MNA circuit models use this
% ROI = [-1 5];  axis3d = [0 5 -1 5];
tfunc_tol = 0.01;
force_deflation = false;

tfunc_index = [2 2];   % which component of a MIMO model transfer function to plot
use_available_URM_tfunc = false;
tfunc_plotstyle = 'linear'; % or 'log'
output_eig_table = true;
plot_tfunc_surface = true; plotres = [400 400]; 
plot_implicit_tfunc_poles = false;
spyH_size = 15;
outputlatex = false;
pause_after_cycle = false;


% Model setup
model_name = sprintf('%s',data_filename);
[A E B C FRdomain URM_FR] = realization(data_filename);  % get model data from file
N = size(A,1);

if use_available_URM_tfunc
	if isempty(FRdomain)
		[s FRdomain] = tfdomain(200,ROI);
		URM_FR = URM_freq_response(data_filename,[],[],FRdomain);
	else
		s = 2i*pi*FRdomain;
		ROI = [min(log10(FRdomain))  max(log10(FRdomain))];
		axis3d = [0  ROI(2) ROI];
	end
else
	[s FRdomain] = tfdomain(200,ROI);
end

if isempty(S0)
	S0 = 2*pi*10^(ROI(1) + (ROI(2)-ROI(1))/2);
end

VROM =  zeros(N,0);
Y = zeros(N,0);
T = [];
flops = 0;

for j = 1:length(S0)
	s0 = S0(j);
	ell(j) = size(Y,2);
	[multH R] = make_SI_op(A,E,B,s0);
	m = size(R,2) + ell(j);
	fprintf('\ncycle %d expanding at s_0 = %s,  ell + m = %d + %d \n',j,s0string(s0),ell(j),size(R,2));
	
	if ell(j) > 0
		result = band_Arnoldi(ell(j),[Y R],multH);  % Process the initial band
		if force_deflation
			result.Vh_defl(:,1:ell(j)) = 0;  % force deflation of thick-started Ritz-vectors
		end
	else
		result = band_Arnoldi(m,R,multH);  % Process the initial band
	end
	
	
	n0 = size(result.V,2);
	result = band_Arnoldi(J(j)+ell(j),[Y R],multH,[],n0,result);
	n = size(result.V,2);

% 	output structure matrix if that option is enabled
	if spyH_size > 0
		k = min(spyH_size,n);
		fprintf('%d x %d map of H_%d,   ell+m= %d+%d,  forced deflation: %d\n',...
			k,k,j,ell(j),size(R,2),force_deflation);
		[dummy Hstring Hlatex] = spyH(result.H(1:k,1:k),norm(result.H));
		
		[dummy rhostring rholatex] = spyH(result.rho(1:m,:),norm(result.H));		
		if outputlatex
			disp([Hlatex rholatex])
		else
			disp([Hstring rhostring])
		end
	end
	
	

% 	count flops
	if isreal(s0)
		flops = flops + result.flops;
	else
		flops = flops + 4 * result.flops;
	end
	
	
% 	Determine Ritz-values for this cycle
	[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
	rho2 = result.rho(:,(ell(j)+1):m);
	cV = C'*result.V;
	[mu wt] = ROM_poles(cV, W, rho2, lambda, s0, FRdomain);
	
	
	converged = rr' < ctol;
	ninf = d2S(mu,s) < 1e10;
	in_ROI = imag(mu) > ROI(1) & imag(mu) < ROI(2);
	plottable = ninf;
	
	
% 	Set criteria for Ritz vectors to keep
	rel_wt = wt / sum(wt);
	keep = (rr' < keep_tol | rel_wt > rel_wt_tol);
	
	if output_eig_table
% 		Output table of 10 most converged eigs
		[rrsort idx] = sort(rr);
		if outputlatex
			fprintf('&( $\\Re(\\mu)$, $\\Im(\\mu)$ ) &\\texttt{rr} & wt & keep \\\\ \n \\hline \n')
			for i = 1:min(10,length(mu))
				fprintf('%d & (\\texttt{% 5.3e}, \\texttt{% 5.3e}) & %g & %g & %d \\\\\n',...
					i, real(mu(idx(i))),imag(mu(idx(i))),rr(idx(i)), wt(idx(i)), keep(idx(i)) );
			end
			fprintf('\\hline\n');
		else
			for i = 1:min(10,length(mu))
				fprintf('%2d  (% 5.3e,% 5.3e )    %4.2e     %4.2e     %1d\n',...
					i, real(mu(idx(i))),imag(mu(idx(i))),rr(idx(i)), wt(idx(i)), keep(idx(i)) );
			end
			fprintf('\n');		
		end
	end
	
	ntot = size(VROM,2)+n-ell(j);
	fprintf('... ROM: %d, n_tot: %d,   converged: %d,   keep: %d  weight: %g\n',...
		n,ntot,nnz(converged), nnz(keep), sum(wt) );
	
	fprintf('...updating thick-restart basis...')
	Znew = result.V*W(:,keep);
	[Y T] = add2basis(Y, T, Znew,dtol);
	fprintf('dim Y = %d\n',size(Y,2));
	
	if plot_implicit_tfunc_poles
		% Plot finite implicit ROM poles
		figure('name',sprintf('poles %d',j))
		plot_ROMpoles(mu(plottable), s0, rr(plottable), wt(plottable), ROI,keep(plottable));
		if j < length(S0)
			hold on
			plot(real(S0(j+1)), imag(S0(j+1)), 'kp');
			hold off
		end
		h1 = gca;
		putaxis;
		title(sprintf('n_{ROM}: %d, weight: %g,   n_{tot}: %d, cvrgd: %d, to keep: %d',...
			n, sum(wt), ntot, nnz(converged), nnz(keep)));
		colorbar;
	end
	
	
	
	VROM = [VROM result.V(:,ell(j)+1:n)];
	j = j+1;
	
	if pause_after_cycle
		input('');
	end
end % of j-th cycle of outer loop


%%%%%%%%%    Post basis-building analysis  %%%%%%%%%%%%%%
[V VROM_split] = make_basis_real(VROM);
nreal = size(V,2);
ROM_exp_FR = transfer_function(V,A,E,C,B,s);


% Determine ROM error or estimate
if use_available_URM_tfunc
	tf_err_strng = sprintf('%g',tf_err);
else
% 	 Output ROM stats to the console
	ROM_exp_FR_previous_iter = transfer_function(V(:,1:(nreal-10)),A,E,C,B,s);
	err_est = tfunc_err(ROM_exp_FR,ROM_exp_FR_previous_iter);
	tf_err_strng = sprintf('%g (est)',err_est);
end

% plot ROM frequency response
for i = 1:size(tfunc_index,1)
	Hi = tfunc_index(i,1); Hj = tfunc_index(i,2);
	figure('name',sprintf('FR (%d,%d)',Hi,Hj));
	ROM_tfunc = abs(squeeze(ROM_exp_FR(Hi,Hj,:)));
	
	if use_available_URM_tfunc
		URM_tfunc = abs(squeeze(URM_FR(Hi,Hj,:)));
		if strcmp(tfunc_plotstyle, 'linear')
			semilogy(FRdomain,URM_tfunc,'r--',FRdomain,ROM_tfunc)
		else
			loglog(FRdomain,URM_tfunc,'r--',FRdomain,ROM_tfunc)
		end
		legend('URM','ROM')
	else
		if strcmp(tfunc_plotstyle, 'linear')
			semilogy(FRdomain,ROM_tfunc)
		else
			loglog(FRdomain,ROM_tfunc)
		end
	end
	title(sprintf('%s: %s',model_name,tf_err_strng));
end

[Usplit singVsplit] = svd(VROM_split,'econ');
c2 = cond(VROM_split);
er2 = rank(VROM_split);
singVsplit = diag(singVsplit);





%  Output ROM stats to the console
fprintf('iterations: %d, ROM size: %d, eff-rank: %g, rel-error: %s,  flops: %d + %dM\n\n',...
	size(VROM,2),nreal,er2,tf_err_strng,flops,length(J));


if use_available_URM_tfunc  % more detailed error analysis
	Compute tfunc-errors vs singular value
	h = waitbar(0,'constructing err vs \sigma...');
	abort = false;
	for n = 1:nreal
		ROM_exp_FR = transfer_function(Usplit(:,1:n),A,E,C,B,s);
		tf_err_sing(n) = tfunc_err(URM_FR,ROM_exp_FR);
		if ~ishandle(h)
			abort = true;
			break;
		else
			waitbar(n/nreal,h);
		end
	end
	if ~abort
		delete(h);
		figure('name','svd2');
		semilogy(1:nreal,singVsplit,'.', 1:nreal,tf_err_sing,'+');
		title(sprintf('cond(V_{ROM}^*): %g,  LI: %g',c2,er2));
		legend('\sigma','rel-err')
	end
	
	
	
% 	Compute tfunc-errors for truncated sub-bases of V from 1 to n
	h = waitbar(0,'constructing err vs n...');
	abort = false;
	for n = 1:nreal
		ROM_exp_FR = transfer_function(V(:,1:n),A,E,C,B,s);
		[tf_err(n) tf_full(n,:)] = tfunc_err(URM_FR,ROM_exp_FR);
		if ~ishandle(h)
			abort = true;
			break;
		else
			waitbar(n/nreal,h);
		end
	end
	
	if ~abort
		delete(h);
		
% 		 plot full-error (fire plot)
		figure('name','fire')
		hp = plotfullerr(tf_full,tfunc_tol);
		title(sprintf('%s: explicit ROM full err',model_name));
		
% 		plot frequency response norm-error
		figure('name','tf_err');
		semilogy(1:nreal,tf_err(1:n),'k.');
		line([1 nreal],[tfunc_tol tfunc_tol],'LineStyle',':');
		title(sprintf('%s: explicit ROM',model_name));
	end
	
end

%  Explicitly project system realization (A,E,B,C) onto V
An = V'*A*V;  En = V'*E*V; Bn = V'*B; Cn = V'*C;

% Plot explicit ROM poles
[Z Mu] = eig(An,En);
mu_exp = diag(Mu);

% compute residues/weights of explicit model
Bh = colnorms((Z\Bn).');
Ch = colnorms(Cn'*Z);
residue = Bh(:) .* Ch(:);
residue(isinf(mu_exp)) = 1;
wt_exp = abs(residue) ./  (1+d2S(mu_exp,FRdomain));
Z = V*Z;
AZ = A*Z;
rr_exp = colnorms(A*Z - E*Z*Mu) ./ colnorms(AZ);
near_S =  d2S(mu_exp,s) < 1e11;

% Plot poles of explicitly projected ROM
figure('name','exp_poles')
plot_ROMpoles(mu_exp(near_S),S0,rr_exp(near_S),wt_exp(near_S),ROI,false(1,nnz(near_S)));
caxis([-3 0])
h2 = gca;
putaxis;

if plot_implicit_tfunc_poles
	linkaxes([h1 h2]);
end

% % try to match poles of the two different ROMs
% for i = 1:length(mu)
% 	[pole_diff idx] = min(abs(mu(i)-mu_exp));
% 	rel_pole_diff(i) = pole_diff / abs(mu(i));
% 	mu2(i) = mu_exp(idx);
% 	rr2(i) = rr_exp(idx);
% end
% figure;
% % disp([abs(mu(:)) rr(:) abs(mu2(:)) rr2(:) rel_pole_diff(:)]);
% loglog(rel_pole_diff,rr(:),'o',rel_pole_diff,rr2(:),'.');
% % line([1 nn],[ctol ctol],'LineStyle',':');
% legend('implicit','explicit');
% title(sprintf('%s ROM, n=%d: FR, err=%g',model_name,size(V,2),tf_err(n)));



%  plot transfer function surface
if plot_tfunc_surface
	if exist('axis3d','var')
		for i = 1:size(tfunc_index,1)
			Hi = tfunc_index(i,1); Hj = tfunc_index(i,2);
			rzn.A = An;
			rzn.E = En;
			rzn.B = Bn(:,Hi);
			rzn.C = Cn(:,Hj);
			abort = tfunc_surf_exp(rzn,S0,plotres,axis3d,ROI);
			if abort
				break;
			else
				title(sprintf('%s',model_name));
				set(gcf,'name',sprintf('surf (%d,%d)',Hi,Hj));
				axis off
				zoom on
			end
		end
	end
end




	function hp = plotfullerr(ferr,tfunc_conv_tol)
		hp = contourf(FRdomain,(1:size(ferr,1)).',log10(ferr),32); % use contourf or pcolor
		set(gca,'YDir','normal');
		shading flat;
		caxis([log10(tfunc_conv_tol) 0]);
		colormap('hot');
		colorbar;
	end


	
end % main

