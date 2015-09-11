function [result dummy] = nrmods(datafilename,S,J,nsamps,snaps)
% This is code to analyze and compare transfer functions, pole convergence, etc
%  of multiple ROMs from a given set of example data.
%  ** This version only computes complex arnoldi iterations and allows for restarts.
%     nmods uses a single run and computes equivalent-real forms
pole_conv_tol = 1e-6;  % convergence tolerance for eigs/poles
tfunc_conv_tol = 1e-5;
sing_tol = sqrt(eps);
nfirst = 2;

dummy = [];

do_tf = true;
do_pol = false;


[C G c b] = inputdata(datafilename);

tfunc_URM = examp_tfunc(datafilename);
[N m] = size(b);
[s frq] = tfdomain();

%% compute Arnoldi iterations
fprintf('performing Rational Arnoldi iterations: \n')
[V H] = RA(A,E,b,S,J);
nmax = size(V,2);

if ~isreal(V)
	V = make_basis_real(V);
end
n = size(V,2);

if ~exist('nsamps','var') || nsamps > nmax-nfirst +1
    nsamps = nmax-nfirst+1;
end

if nsamps < nmax
    nvals = fix(linspace(nfirst,nmax,nsamps));
else
    nvals = nfirst:nmax;
end

Lnv = numel(nvals);
if Lnv == 1
    snaps = nmax;
end

numtrials = 2;
errs = zeros(Lnv,1);
n_cvrgd_pols = zeros(Lnv,numtrials);
mean_rr = zeros(Lnv,numtrials);
med_rr = zeros(Lnv,numtrials);
fullerr = ones(200,Lnv,numtrials);
ROMsize = zeros(Lnv,1);

if ~exist('snaps','var')
    snaps = [];
end
snaps = [snaps(:); nmax+1];
snapidx = 1;


%% Main Loop
% Take convergence samples for each value in nvals  and plot data for each
%  value of n in snaps

for j = 1:Lnv
    n = nvals(j);
%     fprintf(' %d',n);
    if do_cplx
        %% Compute models via complex Krylov-subspace projection
        Vn = orth(V(:,1:n));
%         Vn = V(:,1:n);
        if j == 1 || size(Vn,2) ~= cplx_ROMsize(j-1)
            
            if do_pol
                %% Complex Projection & determine Ritz-poles
                [Cn Gn cn bn mu_pc rr_pc] = project(Vn,C,G,c,b,sing_tol);
                n_cvrgd_pols(j,1) = nnz(rr_pc <= pole_conv_tol);
                mean_rr(j,1) = 10^mean(log10(rr_pc));
                med_rr(j,1) = median(rr_pc);
            else
                %% Complex Projection only
                [Cn Gn cn bn] = project(Vn,C,G,c,b,sing_tol);
            end
            
            ROMsize(j) = length(bn);
            
            if do_tf
                tfunc_pc = transfer_function(Cn,Gn,cn,bn,true);
                [errs(j,1) fullerr(:,j,1)] = tfunc_err(tfunc_pc.',tfunc_URM.');
            end
            
            
        else
            ROMsize(j) = ROMsize(j-1);
            errs(j,1) = errs(j-1,1);
            fullerr(:,j,1) = fullerr(:,j-1,1);
            mean_rr(j,1) = mean_rr(j-1,1);
            med_rr(j,1) = med_rr(j-1,1);
            n_cvrgd_pols(j,1) = n_cvrgd_pols(j-1,1);
        end
    end
    
   
    
    %% Handle Snapshots
    if n >= snaps(snapidx)
        if do_tf
                approx_tfuncs = [tfunc_pc; tfunc_pr];
                %% Plot transfer function(s) at snapshot point
                figure('name',sprintf('tfunc, n=%d',n));
                plot_tfunc(approx_tfuncs,tfunc_URM);
                legend('URM',sprintf('cplx (%d),%3.2g',ROMsize(j),errs(j,1)));
                title(sprintf('%s tfunc, n=%d',datafilename,n));
                
                %% Detailed tfunc error plot (rel error vs S)
                figure('name',sprintf('tfunc err, n=%d',n));
                semilogy(frq,squeeze(fullerr(:,j,:)));
                legend(sprintf('cplx (%d),%3.2g',ROMsize(j),errs(j,1)));
                title(sprintf('%s relerr, n=%d',datafilename,n));
            end
        end
        
        if do_pol
            %% Plot relative residuals for proj ROM Ritz values
            figure('name',sprintf('proj rr, n=%d',n));
            semilogy(sort(rr_pc),'o');
           
            if any(n_cvrgd_pols(j,:))
                line([1 length(rr_pc)],[pole_conv_tol pole_conv_tol],'linestyle','--');
            end
            hold off;
            legend('rr_{cplx proj}','rr_{eqreal proj}','conv tol');
            title(sprintf('relative residuals, n=%d',n));
            
            %% Plot poles in the complex plane
            figure('name',sprintf('poles, n=%d',n));
            plot_pols;
        end
        
        drawnow;
        snapidx = snapidx + 1;
    end
end

fprintf('\n');

result.errs = errs;
result.n_cvrgd_pols = n_cvrgd_pols;
result.mean_rr = mean_rr;
result.med_rr = med_rr;
result.fullerr = fullerr;
result.ROMsize = ROMsize;
result.s0vals = s0vals;

if Lnv > 1 && nargout ~= 1
    % plot model size vs number of iterations
%     figure('name','size vs n');
%             plot(nvals,ROMsize);
%     title('proj ROM size vs n (iters)');
    
    if do_tf
        %% plot full transfer function error vs n contours
            plotfullerr(fullerr(:,:,1),'cplx proj');
       
        
        %% plot norm transfer function err vs n
        figure('name','tfunc errs');
        errs(errs>1) = 1;
        h = semilogy(nvals,errs);
        setlegend(h);
        title('tfunc errors')
    end
    
    if do_pol
        %% plot # of converged poles vs n
        figure('name','cvrgd poles');
        h = plot(nvals,n_cvrgd_pols);
        
        % plot mean relative residual vs n
        figure('name','mean rr');
        h = semilogy(nvals,mean_rr);
        setlegend(h);
        
        % plot median relative residual vs n
        figure('name','med rr')
        h = semilogy(nvals,med_rr);
        setlegend(h);
    end
    
    if do_pol && do_tf
        %  This is some analysis to use to possibly link current status of
        % relative residuals with quality of the ROM computed thus far,
        % allowing us to avoid computing intermediate transfer functions.
        %  Better yet, links between convergence of the eqreal model
        %    and projected model would be nice.
        
        
        %% Superimposed plot of # converged poles & tfunc error vs n
        %         figure('name','proj err & #cvrgd')
        %         [ax h1 h2] = plotyy(nvals,errs(:,1),nvals,n_cvrgd_pols(:,1),@semilogy,@plot);
        %         set(h1,'Marker','o')
        %         set(h2,'Marker','+')
        %         legend('err','# cvrgd')
        
        %% Mean rr vs tfunc err
        %         figure('name','mean rr vs tf err');
        %         h = loglog(mean_rr,errs,'.');
        %         xlabel('mean rr');
        %         ylabel('err')
        
        
        %% # converged poles vs tfunc err
        %         figure('name','tf err vs #cvrgd');
        %         h = semilogy(n_cvrgd_pols,errs,'.');
        %         xlabel('# cvrgd');
        %         ylabel('err')
    end
end

%% ------------------------------------------------------------------------

    function hp = plotfullerr(ferr,txt)
        figure('name',sprintf('%s fullerr',txt));
        hp = contourf(frq,nvals,log10(ferr).',32); % use contourf or pcolor
        title(sprintf('%s err vs f vs n',datafilename));
        set(gca,'YDir','normal');
        shading flat;
        caxis([log10(tfunc_conv_tol) 0]);
        colormap('hot')
        colorbar;
        title(txt)
    end

    function h = plot_pols()
        % disregard infinite poles
        mu_inf = isinf(mu_pc);
        mu_pc(mu_inf) = [];
        rr_pc(mu_inf) = [];
        rmupc = lscale(real(mu_pc));
        imupc = imag(mu_pc);
        
        mu_inf = isinf(mu_pr);
        mu_pr(mu_inf) = [];
        rr_pr(mu_inf) = [];
        rmupr = lscale(real(mu_pr));
        imupr = imag(mu_pr);
        
        
        rs0 = lscale(real(s0vals));
        is0 = imag(s0vals);
        
        s = getS;
        s1 = imag(s(1));
        s2 = imag(s(end));
        
        hold on
        h1 = scatter(rmupc,imupc,[],log10(rr_pc));
        h2 = scatter(rmupr,imupr,[],log10(rr_pr));
        h = [h1 h2];
        [lh oh] = setlegend(h);
        ph=findobj(oh,'type','patch');
        set(ph,'MarkerEdgeColor','k'); % this keeps the legend entries visible
        
        scatter(rs0,is0,225,'kx');
        line([0 0],[s1 s2],'lineWidth',2,'LineStyle','--');
        hold off
        caxis([log10(pole_conv_tol) 0]);
        colormap('autumn');
        colorbar;
	end
end
