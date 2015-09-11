function errs = nmods(datafilename,nmax,nsamps,s0,snaps)
% This is code to analyze and compare transfer functions, pole convergence, etc
%  of multiple ROMs from a given set of example data.

defl_tol = sqrt(eps);  % deflation tolerance for orthogonalization (in Arnoldi and QR)
defl_tol_gs = defl_tol;
pole_conv_tol = 1e-5;
tfunc_conv_tol = 1e-5;

do_tf = false;
do_pol = true;

if ischar(datafilename)
    [C G c b] = inputdata(datafilename);
    tfunc_URM = tfunc_urm(datafilename,C,G,c,b);
    N = length(C);
    [A R0] = makeA(C,G,b,s0);
    [Ar Rr0] = makeAr(C,G,b,s0);
    cr = [c; -1i*c];
else
    % generate random integer test matrices instead of file data
    N(1) = datafilename;
    if numel(datafilename) > 1
        m = datafilename(2);
    else 
        m = 1;
    end
    [A Ar R0 Rr0] = testmat(N,m);
end
    
if nsamps < nmax
    nvals = fix(linspace(2,nmax,nsamps));
else
    nvals = 2:nmax;
end

rho = [norm(R0); zeros(nmax-1,1)];
[s frq] = getS;

tol.defl_flag = 1;
tol.defl_tol = defl_tol;
tol.defl_tol_gs = defl_tol_gs;
tol.normA_flag = 1;
tol.normA = 0;



%% compute Arnoldi iterations
fprintf('performing Arnoldi iterations: ')
% %eqreal Arnoldi
fprintf('eqreal...')
result = paired_band_Arnoldi_mod(nmax,Rr0,Ar,tol);
Hr = result.H;
Vr = result.V;
W = result.W;

Vrcount = result.Vcount;
Wcount = result.Wcount;
vrnext = result.Vh_defl;
nmax_Vr = result.n0;


% standard complex Arnoldi
fprintf('complex...')
result2 = band_Arnoldi(nmax,R0,A,tol);
Hc = result2.H;
Vc = result2.V;
vcnext = result2.Vh_defl;
nmax_Vc = result2.n0;

%%  Comparing PCC Krylov subspace obtained via different routes
% Vrsplit = reshape(Vr,N,[]);
% Vcsplit = reshape([real(Vc); imag(Vc)],N,[]);
% 
% Wr = orth(Vrsplit);
% Wc = orth(Vcsplit);
% [Wr U1] = qr(Vrsplit,0); 
% [Wc U2] = qr(Vcsplit,0); 


%%
numtrials = 3;
Lnv = length(nvals);
errs = zeros(Lnv,numtrials);
n_cvrgd_pols = zeros(Lnv,numtrials);
mean_rr = zeros(Lnv,numtrials);
med_rr = zeros(Lnv,numtrials);
fullerr = ones(200,Lnv,numtrials);

if ~exist('snaps','var')
    snaps = [];
end
snaps = [snaps(:); nmax+1];
snapidx = 1;

lastnVr = 0;
lastnW = 0;
lastnc = 0;

%% Take convergence samples for each value in nvals
%  and plot data for each value of n in snaps
fprintf('\ncomputing stuff: n=');

for j = 1:Lnv
    n = nvals(j);
    fprintf(' %d',nvals(j));
    nr = min(n, nmax_Vr); % # of Vr vectors that span n-th iterative space
    nVr = Vrcount(nr);
    if lastnVr ~= nVr
        lastnVr = nVr;
               
        Vr_n = Vr(:,1:nVr);
        Hr_n = Hr(1:nVr,1:nVr);
        rhor_n = rho(1:nVr);
    
        nW = Wcount(nr);   % # of W vectors that span n-th iterative space
        % no need to re-do projection if PCC-Krylov basis hasn't grown
        if lastnW == nW
            do_proj = false;
            fprintf('\nsaved some comps!\n')
        else
            lastnW = nW;
            do_proj = true;
            splitVr_n = W(:,1:nW);
            sz(j) = size(splitVr_n,2);
        end
        
        if do_tf
            % compute eqreal and proj transfer functions
            tfunc_hess_eqrl = tf_hess(Vr_n,Hr_n,rhor_n,cr,s0); % tfunc via eqreal Hess
            if do_proj
                tfunc_proj = tf_proj(splitVr_n,C,G,c,b); % tfunc via eqreal exp projection
            end
            approx_tfuncs = [tfunc_hess_eqrl; tfunc_proj];
            [errs(j,2:3) fullerr(:,j,2:3)] = tfunc_err(approx_tfuncs.',tfunc_URM.');
        end
        if do_pol
            % determine Ritz-poles via Hr_n
            [mu_r W_r rr_r] = approx_poles_hess(Vr_n,Hr_n,vrnext,cr,rhor_n,s0);
            n_cvrgd_pols(j,2) = nnz(rr_r <= pole_conv_tol);
            mean_rr(j,2) = 10^mean(log10(rr_r));
            med_rr(j,2) = median(rr_r);
            
            % determine Ritz-poles via proj
            if do_proj
                [mu_proj Z_proj rr_proj] = approx_poles_proj(splitVr_n,C,G,c,b,s0,[],'e');
            end
            n_cvrgd_pols(j,3) = nnz(rr_proj <= pole_conv_tol);
            mean_rr(j,3) = 10^mean(log10(rr_proj));
            med_rr(j,3) = median(rr_proj);
        end
    else
        % no point in recomputing transfer functions and errors if nr is the
        % same as before
        sz(j) = sz(j-1);
        errs(j,2:3) = errs(j-1,2:3);
        fullerr(:,j,2:3) = fullerr(:,j-1,2:3);
        n_cvrgd_pols(j,2:3) = n_cvrgd_pols(j-1,2:3);
        mean_rr(j,2:3) = mean_rr(j-1,2:3);
        med_rr(j,2:3) = med_rr(j-1,2:3);
    end
    
    nc = min(n, nmax_Vc);
    if lastnc ~= nc
        lastnc = nc;
        Vc_n = Vc(:,1:nc);
        Hc_n = Hc(1:nc,1:nc);
        rhoc_n = rho(1:nc);
        
        if do_tf
            tfunc_hess_cplx = tf_hess(Vc_n,Hc_n,rhoc_n,c,s0);  % tfunc via complex Hess
            [errs(j,1) fullerr(:,j,1)] = tfunc_err(tfunc_hess_cplx.',tfunc_URM.');
            approx_tfuncs = [tfunc_hess_cplx; approx_tfuncs];
        end
        if do_pol
            % determine Ritz-poles via Hc_n
            [mu_c W_c rr_c] = approx_poles_hess(Vc_n,Hc_n,vcnext,c,rhoc_n,s0);
            n_cvrgd_pols(j,1) = nnz(rr_c <= pole_conv_tol);
            mean_rr(j,1) = 10^mean(log10(rr_c));
            med_rr(j,1) = median(rr_c);
        end
    else
        % no point in re-computing poles, transfer functions,and errors if
        % nc is the same as before
        errs(j,1) = errs(j-1,1);
        fullerr(:,j,1) = fullerr(:,j-1,1);
        n_cvrgd_pols(j,1) = n_cvrgd_pols(j-1,1);
        mean_rr(j,1) = mean_rr(j-1,1);
        med_rr(j,1) = med_rr(j-1,1);
    end
    
    
    %% Analysis at individual n values
    if n >= snaps(snapidx)
        if do_tf %&& false
            %% Plot transfer function at snapshot point
            figure('name',sprintf('tfunc, n=%d',n));
            plot_tfunc(approx_tfuncs,tfunc_URM);
            title(sprintf('%s tfunc, n=%d, s0=%3.2g+%3.2gi',...
                datafilename,n,real(s0),imag(s0)));
            legend('URM',sprintf('hess cplx,%3.2g',errs(j,1)),...
                sprintf('hess eqrl, %3.2g',errs(j,2)), sprintf('proj, %3.2g',errs(j,3))  );
            
            %% Detailed tfunc error plot (rel error vs S)
            figure('name',sprintf('tferr, n=%d',n));
            semilogy(frq,squeeze(fullerr(:,j,:)));
            title(sprintf('%s relerr, n=%d, s0=%3.2g+%3.2gi',datafilename,n,real(s0),imag(s0)));
            legend(sprintf('cplx (hess) ,%3.2g',errs(j,1)),sprintf('eqrl (hess), %3.2g',errs(j,2)),...
                sprintf('proj, %3.2g',errs(j,3))  );
            
            %% Plot tfunc err vs distance to s0
            ds0 = abs(getS - s0);
            figure('name',sprintf('tf err vs dist, n=%d',n));
            h = semilogy(ds0.',squeeze(fullerr(:,j,:)),'.');
            setlegend(h);
        end
        
        if do_pol
            %% Plot relative residuals for (proj) ROM Ritz values
            figure('name',sprintf('proj rr, n=%d',n));
            semilogy(sort(rr_proj),'.');
            if n_cvrgd_pols(j,3)
                hold on
                line([1 length(rr_proj)],[pole_conv_tol pole_conv_tol],'linestyle','--');
                hold off
                legend('rr_{proj}','conv tol');
            end
            title(sprintf('rr_{proj}, sz=%d',sz(j)));
            
            
            %% Plot relative residuals of poles vs their distances from s0
            figure('name',sprintf('rr vs dist, n=%d',n));
            h = loglog(abs(mu_c-s0),rr_c,abs(mu_r-s0),rr_r,abs(mu_proj-s0),rr_proj);
            set(h,'LineStyle','none');
            title('rr_{\mu} vs |\mu- s_0|');
            setlegend(h);
            
            %% Plot distribution of relative residuals (histogram)
            figure('name',sprintf('rr distrib, n=%d',n));
            bins = linspace(-10,0,11)';
            nc = hist(log10(rr_c),bins);
            nr = hist(log10(rr_r),bins);
            np = hist(log10(rr_proj),bins);
            bar(bins,[nc; nr; np].','group');
            legend('cplx (Hess)','eqrl (Hess)','proj')
            
            %% Plot poles in the complex plane
            figure('name',sprintf('poles, n=%d',n));
            plot_pols;
%             axis(good_pole_domain);
        end
        
        drawnow;
        snapidx = snapidx + 1;
    end
end
fprintf('\n')

if Lnv > 1
    %% plot projected ROM size vs n
        figure('name','size vs n');
        plot(nvals,sz,1:nmax_Vr,Wcount);
        title('proj ROM size vs n (iters)')
    
    if do_tf
        %% plot full transfer function error vs n contours
        plotfullerr(fullerr(:,:,2),'eqreal');
        plotfullerr(fullerr(:,:,3),'proj');
        
        %% plot norm transfer function err vs n
        figure('name','errs');
        h = semilogy(nvals,errs);
        tt = ' tfunc err';
        title(sprintf('%s %s, s0=%1.1g+i%1.1g',datafilename,tt,[real(s0) imag(s0)]));
        setlegend(h);
    end
    
    if do_pol
        %% plot # of converged poles vs n
        figure('name','cvrgd poles');
        h = plot(nvals,n_cvrgd_pols);
        setlegend(h);
        tt = ' #cvrgd';
        title(sprintf('%s %s, s0=%1.1g+i%1.1g',datafilename,tt,[real(s0) imag(s0)]));
        
        %% plot mean relative residual vs n
        figure('name','mean rr');
        h = semilogy(nvals,mean_rr);
        title(sprintf('%s, s0=%1.1g+i%1.1g',datafilename,[real(s0) imag(s0)]));
        setlegend(h);
        
        %% plot median relative residual vs n
        figure('name','med rr')
        h = semilogy(nvals,med_rr);
        title(sprintf('%s, s0=%1.1g+i%1.1g',datafilename,[real(s0) imag(s0)]));
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
%         [ax h1 h2] = plotyy(nvals,errs(:,3),nvals,n_cvrgd_pols(:,3),@semilogy,@plot);
%         set(h1,'Marker','o')
%         set(h2,'Marker','+')
%         legend('err','# cvrgd')
        
        %% Mean rr vs tfunc err
%         figure('name','mean rr vs tf err');
%         h = loglog(mean_rr,errs,'.');
%         title(sprintf('%s, s0=%1.1g+i%1.1g',datafilename,[real(s0) imag(s0)]));
%         xlabel('mean rr');
%         ylabel('err')
%         setlegend(h);
        
         %% # converged poles vs tfunc err
%         figure('name','tf err vs #cvrgd');
%         h = semilogy(n_cvrgd_pols,errs,'.');
%         title(sprintf('%s, s0=%1.1g+i%1.1g',datafilename,[real(s0) imag(s0)]));
%         xlabel('# cvrgd');
%         ylabel('err')
%         setlegend(h);
        
        %% eqreal mean rr vs proj tfunc err
%         figure('name','proj tf err vs eqrl mean rr');
%         h = loglog(mean_rr(:,2),errs(:,3),'.');
%         title(sprintf('%s, s0=%1.1g+i%1.1g',datafilename,[real(s0) imag(s0)]));
%         xlabel('eqrl mean rr');
%         ylabel('proj err')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function er = formerr(Z)
%         Zt = Z(1:N,:);
%         Zb = Z(N+1:end,:);
%         er = (abs(1/2 + sum([-imag(Zt); real(Zt)] .* [real(Zb); imag(Zb)])) ).';
%     end

    function hp = plotfullerr(ferr,txt)
        figure('name',sprintf('%s fullerr',txt));
        hp = contourf(frq,nvals,log10(ferr).',32); % use contourf or pcolor
        title(sprintf('%s err vs f vs n, s0=%1.2g+i%1.2g',datafilename,[real(s0) imag(s0)]));
        set(gca,'YDir','normal');
        shading flat;
        caxis([log10(tfunc_conv_tol) 0]);
%         colorbar;
        colormap('hot')
    end

    function h = plot_pols
        % disregard infinite poles
        mu_inf = isinf(mu_c);
        mu_c(mu_inf) = [];
        rr_c(mu_inf) = [];
        mu_inf = isinf(mu_r);
        mu_r(mu_inf) = [];
        rr_r(mu_inf) = [];
        mu_inf = isinf(mu_proj);
        mu_proj(mu_inf) = [];
        rr_proj(mu_inf) = [];
        
        rmuc = lscale(real(mu_c));
        imuc = lscale(imag(mu_c));
        rmur = lscale(real(mu_r));
        imur = lscale(imag(mu_r));
        rmup = lscale(real(mu_proj));
        imup = lscale(imag(mu_proj));
        rs0 = lscale(real(s0));
        is0 = lscale(imag(s0));
        s = getS;
        s1 = lscale(imag(s(1))); 
        s2 = lscale(imag(s(end)));
        
        hold on
        h1 = scatter(rmuc,imuc,[],log10(rr_c));
        h2 = scatter(rmur,imur,[],log10(rr_r));
        h3 = scatter(rmup,imup,[],log10(rr_proj));
        h = [h1;h2;h3];
        [lh oh] = setlegend(h);
        ph=findobj(oh,'type','patch');
        set(ph,'MarkerEdgeColor','k'); % this keeps the legend entries visible

        scatter(rs0,is0,225,'kx');
        line([0 0],[s1 s2],'lineWidth',2,'LineStyle','--');
        hold off
        caxis([log10(pole_conv_tol) 0]);
        colormap('hot');
        colorbar;
        
    end

    function [lh oh] = setlegend(h)
        set(h,{'Marker'},{'.','+','o'}');
        [lh oh] = legend(h,'cplx','eqrl','proj');
        
    end
end