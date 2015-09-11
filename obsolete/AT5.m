function [V rel_err] = AT5(datafilename,nvals,s0vals)
% [V rel_err] = AT5(datafilename,nvals,s0vals)
%
% nvals = # iterations for run(s)
% s0vals = (optional) s0 points for each run
%
% plots of ROM poles and transfer function(s) are output if AT5a is called
% with an output argument
%
% example:  AT5('ex308s_1',[50 10 10],(1+1.5i)*pi*1e10)
%
%  this does three runs of thick-restart band_Arnoldi of 50, 10 and 10
%  iterations, starting with a specified initial s0 point and letting
%  the adaptive routine choose the other two expansion points.
%
% example: AT5('ex308s_1',[50 inf 10],[1+1.5i, 1i, 3+2i]*pi*1e10)
%
%  Three runs of thick restart band_arnoldi with the first and third s0
%   predetermined and the second one to be chosen by the adaptive
%   algorithm.

% close all


%% Internally defined parameters
defl_tol = sqrt(eps);  % deflation tolerance for band_Arnoldi
convergence_tol = 1e-5; % poles (eigs) with this relative residual are considered converged
proj_convergence_tol = 6e-3;
conj_tol = 0.1;
wt_good = 0.98;  % for selecting a new s0, consider poles whose combined weight is wt_good
s0init = (0.2+5i)*2*pi*1e9; % default initial s0 (should be generally located)
cconj = false; % if true, combine complex conjugate poles and tf_term pairs
lp = true;    % if true, display poles on log-scaled real axis

% good_pole_domain = [-5e9 5e7 -1.5e10 3*pi*1e10]; % specify axis limits for plotting poles
good_pole_domain = [];

%% Setup

% create vector of s0 values.  Those to be autoselected are indicated by inf.
rmax = size(nvals,2);
S0 = inf(1,rmax);
if exist('s0vals','var') && ~isempty(s0vals)
    S0(1:length(s0vals)) = s0vals;
else
    S0(1) = s0init;
end
s0 = S0(1);

% Input data (matrices C,G,b) from data file; note H(s) = c'*(C-s*G)*b
[C G c b] = inputdata(datafilename);
cr = [c; -1i*c];
[Ar Rr0] = makeAr(C,G,b,s0);
[A R0] = makeA(C,G,b,s0);

% Ar is a function that computes Ar(v) = Ar*v, where Ar is the equivalent
% real formulation of A = (C-s0*G)\b.

[N m] = size(b);

Nr = 2*N;
n = nvals;
n_tot = sum(nvals);

Vr = zeros(Nr,n_tot*m);
Yr = zeros(Nr,n_tot*m);
Y = [];

yk = 0;
r = 0;
k_tot = 0;
%% Thick Restart band Arnoldi iterations
while r < rmax
    % generate n Arnoldi vectors,
    if r == 0
        %  First time there is no data from previous runs
        tol.defl_flag = 1;
        tol.defl_tol = defl_tol;
        tol.normA_flag = 1;
        tol.normA = 0;
        n0 = 0;
        result = band_Arnoldi(n(r+1),Rr0,Ar,tol);
        
        if result.exh_flag > 0
            error('Cannot continue! abort after run %d',r);
        end
    else
        % Subsequent runs avoid re-discovering information in
        % the subspace Y, which we already know is invariant.
        % (all previous y's will be orthogonalized against on this run)
        
        clear result0 tol Iv
        % the way we re-start Arnoldi depends on whether s0 changed
        if s0_changed
            % If the new s0 is significantly different from the last one,
            %  compute new matrix and start vector.
            [Ar Rr0] = makeAr(C,G,b,s0);
            [N m] = size(b);
            
            result0.Vh_defl = Rr0;
            result0.m = m;
            result0.mc = m;
            Iv.ph = 1:m;
            Iv.I = [];
            Iv.pd = [];
            Iv.nd = 0;
            result0.Iv = Iv;
            
            tol.defl_flag = 1;
            tol.defl_tol = defl_tol;
            tol.normA_flag = 1;
            tol.normA = 0;
            result0.tol = tol;
        else
            % otherwise, continue where the last run left off.
            result0.Vh_defl = result.Vh_defl;
            result0.m = result.m;
            result0.mc = result.mc;
            result0.Iv = result.Iv;
            tol = result.tol;
            result0.tol = tol;
        end
        
        % basis vectors for invariant space from last run will be treated
        %  as previous vectors on this run. we won't bother computing
        %  the corresponding H (Hessenberg) matrix entries.
        result0.V = Yr(:,1:yk);
        n0 = yk;     % iterations will start with (yk+1)-th vector
        result0.n0 = n0;
        result0.H = [];
        result0.rho = []; % rho is created new, so this isn't really necessary
        result0.exh_flag = 0;
        
        result = band_Arnoldi(n0+n(r+1),Rr0,Ar,tol,n0,result0);
    end
    
    H = result.H;
    vnext = result.Vh_defl;
    Vtmp = result.V;
    V = Vtmp(:,n0+1:result.n0);
    %  V is the set of equivalent real vectors obtained from the
    %  previous run of Arnoldi.
    
    vm = size(V,2);
    
    % combine them with the vectors in storage
    Vr(:,k_tot+1:k_tot+vm) = V;
    k_tot = k_tot+vm;
    
    r = r+1;
    if r < rmax
        % Extract a spanning set of Y of (nearly) invariant subspace of
        %  V and possibly obtain a better expansion point based on
        %  approximate model so far
        if isinf(S0(r+1))
            [Y next_s0] = deflate(V,H,vnext,convergence_tol);
        else
            Y = deflate(V,H,vnext,convergence_tol);
            next_s0 = S0(r+1);
        end
        s0_changed = chord(s0,next_s0) > eps;
        s0 = next_s0;
        S0(r+1) = s0;
        
        l = size(Y,2);
        if l > 0
            % add Y to the (nearly) invariant subspace already discovered
            %  and orthogonalize it.
            Yr(:,yk+1:yk+l) = Y;
            yk = yk+l;
        end
    end
end

if k_tot < size(Vr,2)
    Vr(:,k_tot+1:end) = [];
end
if yk < size(Yr,2)
    Yr(:,yk+1:end) = [];
end


V = orth(reshape(Vr,N,[]));

%% Analysis of ROM. plots etc.

if nargout > 1
    %     % check passivity
    %     Cp = V'*C*V;
    %     Gp = V'*G*V;
    %     if ~isposreal2(Cp,Gp)
    %         fprintf('ROM t-func is not pos real\n');
    %     else
    %         fprintf('ROM t-func is pos real\n');
    %     end
    
    
    
    %%  evaluate ROM transfer function against URM
    tfunc_unreduced = abs(tfunc_urm(datafilename,C,G,c,b));
    tfunc = abs(tf_proj(V,C,G,c,b));
    rel_err = norm(tfunc_unreduced(:)-tfunc(:),inf)/ norm(tfunc_unreduced(:),inf);
    
    %% plot the poles of this ROM
    [mu Z rr tf_terms] = approx_poles_proj(V,C,G,c,b,s0,cconj);
    plot_poles(mu, tf_terms_weight(tf_terms),rr,S0,ax,lp);
    
    
    %% plot transfer function(s)
    n = size(b,2); p = size(c,2);
    siso = n == 1 & p == 1;  % is this a SISO model?
    if siso
        figure
        plot_tfunc(tfunc,tfunc_unreduced,'reduced','unreduced');
        title(sprintf('n = %d,\t rel err = %g',size(V,2),rel_err));
    else
        h=figure;
        aviobj = avifile('tfunc','quality',90,'fps',5);
        [s1 s2 s3] = size(tfunc);
        for i = 1:s1
            for j = i:s2
                %                 figure
                plot_tfunc(tfunc(i,j,:),tfunc_unreduced(i,j,:),'reduced','unreduced');
                aviobj = addframe(aviobj,getframe(h));
                %                 title(sprintf('(%d,%d)',i,j));
            end
        end
        close(h)
        aviobj = close(aviobj);
        %         title(gcf,sprintf('n = %d',size(V,2)));
    end
end

%% -------------------------------------------------------------
% Deflates the set of vectors V, which is provided in equivalent-real
% form. The invariant space Y is returned in equivalent-real form.
    function [Y better_s0] = deflate(V,H,vnext,convergence_tol)
        lct = convergence_tol;
        Vsplit = orth(reshape(V,N,[]));
        
        [mu_proj Z_proj rr_proj tft_proj] = approx_poles_proj(Vsplit,C,G,c,b,s0,false,'qz');
        cvrgd_proj = rr_proj < proj_convergence_tol;
        n_cvrgd_proj = nnz(cvrgd_proj);
        
        % plot rel-residual err of pole approximations via proj
        figure('name','proj rr');
        semilogy(sort(rr_proj),'.');
        if n_cvrgd_proj
            line([1 length(rr_proj)],[proj_convergence_tol proj_convergence_tol],'linestyle','--');
            legend('rr_{proj}','conv tol');
        end
        title('relresid err of pole approx via projection');
        
        
        
        % determine poles, (Arnoldi relation) rel-residual error via eigs of Hn
        [mu W rr L tft_hess] = approx_poles_hess(V,H,vnext,cr,rho,s0);
        cvrgd = rr <= convergence_tol;
        %
        
        
        % compute ROM and URM transfer functions 
        tfunc_proj_tft = sum(tft_proj);
        tfunc_proj = tf_proj(Vsplit,C,G,c,b);
        tfunc_URM = tfunc_urm(datafilename,C,G,c,b);
        tfunc_hess = tf_hess(V,H,rho,cr,s0);
        tfunc_hess_tft = sum(tft_hess);
        %
        approx_tfuncs = [tfunc_proj; tfunc_hess];
        figure;
        errs = plot_tfunc(approx_tfuncs,tfunc_URM);
        title(sprintf('%s, n=%d, s0=%g',datafilename,n(r),s0));
        legend('URM',sprintf('proj,%3.2g',errs(1)),sprintf('Hess, %3.2g',errs(2)) );
        %
        
        
        
        approx_tfuncs = [tfunc_hess; tfunc_hess_tft ];
        figure;
        errs = plot_tfunc(approx_tfuncs,tfunc_URM);
        title(sprintf('%s, n=%d, s0=%g',datafilename,n(r),s0));
        legend('URM',sprintf('hess,%3.2g',errs(1)),sprintf('hess tft, %3.2g',errs(2)) );
        
        approx_tfuncs = [tfunc_proj; tfunc_proj_tft ];
        figure;
        errs = plot_tfunc(approx_tfuncs,tfunc_URM);
        title(sprintf(' %s, n=%d, s0=%g',datafilename,n(r),s0));
        legend('URM',sprintf('proj, %3.2g',errs(1)),sprintf('proj tft, %3.2g',errs(2)) );
        
        
        Z = V*W;   % Ritz vectors
       
        % explicitly compute rel residual error of Z as eigenspace of Heq
        ArZ = complex(Ar(real(Z)),Ar(imag(Z)));
        ZL = Z*diag(L);
%         rr_exp = max(abs(ArZ - ZL))' ./ max(abs(ArZ))';
        rr_exp = sqrt(sum(abs(ArZ - ZL).^2)).' ./ abs(L);
        
        % construct "structurally correct" Zalt  
        Zt = Z(1:N,:);
        Zb = Z(N+1:end,:);
        Zalt = (Z + 1i*[Zb;-Zt])/2;
        
        % explicitly compute rel residual error of Zalt as eigenspace of Heq
        ArZalt = complex(Ar(real(Zalt)), Ar(imag(Zalt)));
        rralt = sqrt(sum(abs(ArZalt - ZL).^2)).' ./ abs(L);
%         rralt = max(abs(ArZalt - ZL))' ./ max(abs(ArZalt))';
        
        % closeness of ritz vectors to their "structural corrections"
        form_err = (abs(1/2 + sum([-imag(Zt); real(Zt)] .* [real(Zb); imag(Zb)])) ).';
%         relZdiff_inf = max(abs(Z-Zalt))'./max(abs(Z))';
%         relZdiff_2 = (sum(abs(Z-Zalt).^2)).';
        
        % We can determine which conjugate is the correct one now
        correct_conjugate = form_err < conj_tol;     
        
        
        
        % from Z for Heq we determine the analogous space X for H
        X =  (Zt + 1i*Zb ) / sqrt(2);
        
      
        % Explicitly compute rel residual error of X as approx eigenspace of H
        AX = A(X);
        XL = X*diag(L);
        rr_eigH = max(abs(AX - XL))' ./ max(abs(AX))';
        
        % Explicitly compute rel residual error of X as approx eigenspace of (A,E)
        CX = C*X;
        rr_eigAE = max(abs(CX - G*X*diag(mu)))' ./ max(abs(CX))';
        
        % compare Z with Zalt (structurally correct mod)
        [dummy rridx] = sort(rr_exp);  % a sorting of poles in order of rr_exp
        chk = [rr_exp rralt];
        figure; semilogy(chk(rridx,:),'.')
        line([1 length(rr_exp)],[lct lct],'linestyle','--');
        legend('Z','Zalt','conv tol')
        title(sprintf('Z and Z_{alt} as eigenspace of Heq'))
        
        
        
        % Compare X approximations (should be almost identical)
%         chk = [mu rr_eigH rr_eigAE]; 
%         figure; semilogy(chk(rridx,2:end),'.');
%         line([1 length(rr)],[lct lct],'linestyle','--')
%         legend('rr_{eigH}','rr_{eigAE}','conv tol')
%         title('X as approx eigenspace of H, and (A,E) -- relresid err')
       
        
        % *** Compare different measures of pole convergence ***
%         chk = [mu rr_exp  rr  form_err]; 
%         figure('name','err_methds');
%         semilogy(chk(rridx,[2:end]),'.');
%         line([1 length(rr)],[lct lct],'linestyle','--')
% %         line([1 length(rr)],[conj_tol conj_tol],'linestyle','-.')
%         legend('rr_{explicit}','rr_{arnoldi}','form-err',...
%            'conv tol','conj tol');
%         title('H_n Ritz convergence measures');
        
   
        % Plot the poles determined via Arnoldi matrix
        wt_hess = tf_terms_weight(tft_hess);
%         chk = [mu wt_hess relZdiff.'];
        
        figure; plot_poles(mu,wt_hess,form_err,s0,true);
        title('poles from Hess');
        axis(good_pole_domain);
        
        % Plot the poles determined via explicit projection
        wt_proj = tf_terms_weight(tft_proj);
        figure; plot_poles(mu_proj,wt_proj,rr_proj,s0,true);
        title('poles from Proj');
        axis(good_pole_domain);
        
        
        % good is logical index of eigenvalues that are (1) of H and not
        % conj(H), and (2) are sufficiently converged
        good = correct_conjugate & cvrgd;
        
        n_good = nnz(good);
        fprintf('%d converged ritz vecs Hess\n',n_good*2);
        fprintf('%d converged ritz vecs by proj\n',n_cvrgd_proj);
        
        %         [tft wt] = tfterms(V,L,X,c,r0);
        
        
        if false && nargout > 1  % only do this if the calling function wants it
            ncvrgd = correct_conjugate & ~cvrgd;
            
            % select a better s0 based on unconverged poles
            better_s0 = autoselect_s0(mu(ncvrgd),wt(ncvrgd),rr(ncvrgd),wt_good);
            %         figure
            %         plot_poles(mu, wt,rr,better_s0,[],false);
            %         close all
        else
            better_s0 = s0;
        end
        
        
        % plot rresid vs distance from expansion point
%         figure; loglog(rr_exp,abs(mu-s0),'.');
%         xlabel('rr')
%         title('rr vs |mu-s_0|');
%         
%         figure; loglog(rr_proj,abs(mu_proj-s0),'.');
%         title('rr vs |mu_{proj}-s_0|');
%         xlabel('rr')
        
        if n_good > 0
            Zgood  = Z(:,good);
            mugood = mu(good);
            rrgood = rr(good);
            rrgood_exp = rr_exp(good);           
           
            
            
            % plot converged poles (from proj model)
            if  n_cvrgd_proj
                mu_rendundant = [mugood;conj(mugood)];
                figure; scatter(real(mu_rendundant),imag(mu_rendundant),10,'k+');
                hold on;
                scatter(real(mu_proj(cvrgd_proj)),imag(mu_proj(cvrgd_proj)),10,'ro');
                hold off;
                legend('hess','proj')
            end
            
            mumu = [mugood ones(n_good,1)  rrgood_exp;...
                mu_proj(cvrgd_proj) zeros(n_cvrgd_proj,1) rr_proj(cvrgd_proj)];
            chk = sortrows(mumu);
           
        
            Y = [real(Zgood) imag(Zgood)];
        else
            Y = [];
        end
    end

%%
    function s0opt = autoselect_s0(mu,wt,rr,wt_good)
        % Auto selects an expansion point s0 based on given convergence
        % information of poles.
        
        % placement (min distance to poles of greatest weight)
        % eliminate poles with positive real part, negative imag part
        bad = real(mu)> 0 | imag(mu)<0;
        mu(bad) = []; wt(bad) = []; rr(bad) = [];
        
        % eliminate poles that deviate from 1 order of magnitude from most
        % poles. these are 'large' or 'small' meaning they are essentially
        % zero or Inf.
        [small large] = magvals(mu,0.7);
        bad = small | large;
        mu(bad) = []; wt(bad) = []; rr(bad) = [];
        
        wt = wt / sum(wt);
        
        % wt2 = wt ./ rr;    % scale up the weight of more converged poles.
        wt2 = wt.*rr;      % scale down the weight of more converged poles
        wt2 = wt2/sum(wt2);
        % wt2 = wt;            % leave weight independent of convergence
        
        % sort poles by weight
        [wt2 idx] = sort(wt2,'descend');
        mu = mu(idx);
        rr = rr(idx); wt = wt(idx);
        
        cwt2 = cumsum(wt2);
        % chk2 = [mu rr wt wt2 cwt2];
        mwt = find(cwt2 > wt_good,1,'first') ;
        
        s0opt= sum(wt2(1:mwt).*mu(1:mwt));
        % s0opt = imag(s0opt);
        
        % plot_poles(mu(1:mwt), wt2(1:mwt),rr(1:mwt),s0opt,[],false);
        % title('selected')
    end
end % of main function
