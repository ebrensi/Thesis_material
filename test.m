function  [S0 J ell flops] = test(data_filename)

% parameters
dtol = sqrt(eps);
ctol = dtol;
keep_tol = 1e-4;
rel_wt_tol = inf;
ROI = [9 10]; % MNA
Nmax = 20;
stepsize = 4;
tfunc_tol = 0.01;
force_deflation = false;
tfunc_index = [1 1];

URM_known = false;
tfunc_plotstyle = 'log';


% Model setup
model_name = sprintf('%s',data_filename);
[s FRdomain] = tfdomain(200,ROI,'linear');
% URM_FR = URM_freq_response(data_filename,[],[],FRdomain);
[A E B C] = realization(data_filename);
N = size(A,1);

VROM =  zeros(N,0);
Y = zeros(N,0);
T = [];

shift_changed = true;
j = 1;
flops = 0;

S0 = getpoint([-10^ROI(2) 10^ROI(2) 10^(ROI(1)-1) 10^(ROI(2)+1)],ROI);

while j <= length(S0)
	if shift_changed
		s0 = S0(j);
		ell(j) = size(Y,2);
		[multH R] = make_SI_op(A,E,B,s0);
		m = size(R,2) + ell(j);
		fprintf('\nCycle %d expanding at s_0 = %s,  band_size = %d + %d \n',j,s0string(s0),size(R,2),ell(j));
		
		if ell(j) > 0
			result = band_Arnoldi(ell(j),[Y R],multH);  % Process the initial band
			if force_deflation
				result.Vh_defl(:,1:ell(j)) = 0;  % force deflation of thick-started Ritz-vectors
			end
			n0 = size(result.V,2);
			result = band_Arnoldi(m,[Y R],multH,[],n0,result);	
		else
			result = band_Arnoldi(m,R,multH);  % Process the initial band
		end
		
		n = size(result.V,2);
		converged = [];
		nmax = Nmax;
		n_eigs = 1;
		tf_err = 1;
	else
		n_eigs = n_eigs + 1;
		nmax = n + Nmax;
	end
	
	while nnz(converged) < ell(j)+ n_eigs && n <= nmax 
		result = band_Arnoldi(n+stepsize,[Y R],multH,[],n,result);
		if length(result.H) ~= n+stepsize
			n = length(result.H);
			break;
		end
		n = n + stepsize;
		
		
		%% **** Inner-loop anaylsis *****
		[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
		converged = rr < ctol;
		
		% 		cV = C'*result.V;
		% 		rho2 = result.rho(:,(ell(j)+1):m);
		% 		[mu wt] = ROM_poles(cV, W, rho2, lambda, s0, FRdomain);
		% 		ROM_wt(n) = sum(wt);
	end
	
	
	
	% This code executes after every cycle at one point s0(j)
	[lambda, W, rr, Vd, F1, Vh, F2] = band_Arnoldi_Ritz(result);
	rho2 = result.rho(:,(ell(j)+1):m);
	
	
	% Get pole-distribution for the j-th ROM
	cV = C' * result.V;
	[mu wt] = ROM_poles(cV, W, rho2, lambda, s0, FRdomain);
	
	converged = rr' < ctol;
	ninf = d2S(mu,s) < 1e10;
	in_ROI = imag(mu) > 0 & imag(mu) < 1e11;
	plottable = ninf & imag(mu) >= 0;
	
	% Set criteria for Ritz vectors to keep
	rel_wt = wt / sum(wt);
	keep = (rr' < keep_tol | rel_wt > rel_wt_tol);
	% 	chk = [abs(mu) rr' wt rel_wt keep]
	
	ntot = size(VROM,2)+n-ell(j);
	fprintf('... ROM: %d, n_tot: %d,   converged: %d,   keep: %d  weight: %g\n',...
		n,ntot,nnz(converged), nnz(keep), sum(wt) );
	
	
	
	%% Plot finite implicit ROM poles
	plot_ROMpoles(mu(plottable), s0, rr(plottable), wt(plottable), ROI,keep(plottable));
	h1 = gca;
	putaxis;
	title(sprintf('n_{ROM}: %d, weight: %g,   n_{tot}: %d, cvrgd: %d, to keep: %d',...
		n, sum(wt), ntot, nnz(converged), nnz(keep)));
	colorbar;
	
	
	
	
	zoom on;   % use mouse button to zoom in or out
	waitfor(gcf,'CurrentCharacter',13)
	zoom reset
	zoom off
	set(gcf,'currentch','0')
	drawnow;
	
	
	[x, y] = ginput(1);
	lims = xlim;
	
	close;
	if x > lims(2)
		shift_changed = false;
	else
		% if shift changed or we're done then update flop count
		% count flops
		J(j) = n-ell(j);
		if isreal(s0)
			flops = flops + result.flops;
		else
			flops = flops + 4 * result.flops;
		end
		
		% update VROM
		VROM = [VROM result.V(:,ell(j)+1:n)];
		
		if x < lims(1)
			break;
		else
			fprintf('...updating thick-restart basis...')
			Znew = result.V*W(:,keep);
%             Znew = [result.V*W(:,~converged & keep) real(result.V*W(:,converged & keep)) imag(result.V*W(:,converged & keep))];
			[Y T] = add2basis(Y, T, Znew,dtol);
			fprintf('dim Y = %d\n',size(Y,2));
			
			
			S0(j+1) = complex(x,y);
			hold on
			plot(real(S0(j+1)),imag(S0(j+1)),'kp');
			hold off
			figure;
			fprintf('... end cycle %d,  size: %d, keeping: %d \n',j,size(VROM,2),nnz(keep));
			
			j = j+1;
			shift_changed = true;
		end
	end
	
end % of j-th cycle of outer loop


% tfunc_tol = 0.01;  % Now set tfunc_tol back to 0.01 for display purposes
% figure;
% h = subplot(2,1,1);
% semilogy(1:n,tf_exp_err(1:n),'k.')
% line([1 n],[tfunc_tol tfunc_tol],'LineStyle',':');
% title(sprintf('%s: s_0 = %s',model_name,s0string));
% set(h,'XTick',[]);
%
% subplot(2,1,2);
% semilogy(ROM_wt,'b.');
% title('ROM weight');
%


%% Post basis-building analysis
V = make_basis_real(VROM);
nreal = size(V,2);
ROM_exp_FR = transfer_function(V,A,E,C,B,s);

% plot ROM frequency response
figure('name','FR')
ROM_tfunc = abs(squeeze(ROM_exp_FR(tfunc_index(1),tfunc_index(2),:)));
if strcmp(tfunc_plotstyle, 'linear') 
	semilogy(FRdomain,ROM_tfunc)
else 
	loglog(FRdomain,ROM_tfunc)
end
% tf_err = tfunc_err(URM_FR,ROM_exp_FR);
% URM_tfunc = abs(squeeze(URM_FR(tfunc_index(1),tfunc_index(2),:)));
% semilogy(FRdomain,ROM_tfunc, FRdomain, URM_tfunc)
% legend(sprintf('ROM: err %g',tf_err),'URM');

An = V'*A*V;  En = V'*E*V; Bn = V'*B; Cn = V'*C;

%% Plot explicit ROM poles
[AA EE QQ QR Z W] = qz(An,En,'complex');  % note that we have renamed V
alpha = diag(AA); beta = diag(EE);
clear QQ QR
mu_exp = alpha ./ beta;
Mu = diag(mu_exp);


%% compute residues/weights
Bh = colnorms(Bn'*W);
Ch = colnorms(Cn'*Z);
residue = Bh(:) .* Ch(:);
residue(isinf(mu_exp)) = 1;
vec_scaling = diag(W'*En*Z);
wt_exp = abs(residue) ./ abs(vec_scaling) ./  (1+d2S(mu_exp,FRdomain));
Z = V*Z;
rr_exp = colnorms( A*Z - E*Z*Mu ) / norm(A,'fro')';
pos =  d2S(mu_exp,s) < 1e10;
figure
plot_ROMpoles(mu_exp(pos),S0,rr_exp(pos),wt_exp(pos),ROI,zeros(1,length(mu)));
caxis([-2 0])
h2 = gca;
putaxis
linkaxes([h1 h2])
% disp(sortrows([mu_exp(pos) rr_exp(pos)' wt_exp(pos)],2));

tf_err = nan;
% Output ROM stats to the console
fprintf('iterations: %d, ROM size: %d, rel-error: %g,  flops: %d + %dM\n\n',...
	size(VROM,2),nreal,tf_err,flops,length(J));



	
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

