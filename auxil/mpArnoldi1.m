function [Vout Yout S0] = mpArnoldi1(data,nvals,s0vals,convergence_tol)
% [Vout Yout S0] = mpArnoldi1(data,nvals,s0vals,convergence_tol)
%
% nvals = # iterations for run(s)
% s0vals = (optional) s0 points for each run
%
% example:  mpArnoldi1('ex308s_1',[50 10 10],(1+1.5i)*pi*1e10)
%
%  this does three runs of thick-restart band_Arnoldi of 50, 10 and 10
%  iterations, starting with a specified initial s0 point and letting
%  the adaptive routine choose the other two expansion points.
%
% example: mpArnoldi1('ex308s_1',[50 inf 10],[1+1.5i, 1i, 3+2i]*pi*1e10)
%
%  Three runs of thick restart band_arnoldi with the first and third s0
%   predetermined and the second one to be chosen by the adaptive algorithm.

% midx = 0;
% M = getframe;
% close 

%% Internally defined parameters
defl_tol = sqrt(eps);  % deflation tolerance for band_Arnoldi
s0init = pi*1e10; % default initial s0 (should be generally located)
min_keep = 10;

%% Setup
% create vector of s0 values.  Those to be autoselected are indicated by inf.
rmax = size(nvals,2);
S0 = nan(1,rmax);
if exist('s0vals','var') && ~isempty(s0vals)
    S0(1:rmax) = s0vals;
else
    S0(1) = s0init;
end

if ~exist('convergence_tol','var')
    convergence_tol = 1e-5;
end

% Input data (matrices C,G,b) from data file; note H(s) = c'*(C-s*G)*b
[C G c b] = inputdata(data);

[N m] = size(b);
n = nvals;

Vtot = [];

%% Thick Restart band Arnoldi iterations
for r = 1:rmax
    s0 = S0(r);
    
    % generate n Arnoldi vectors,
    if r == 1
        %  First time there is no data from previous runs
        tol.defl_flag = 1;
        tol.defl_tol = defl_tol;
        tol.normA_flag = 1;
        tol.normA = 0;
        n0 = 0;
        [A R0] = makeA(C,G,b,s0);
        result = band_Arnoldi(n(1),R0,A,tol);
        Up = [];
        Y = [];
    else
        % Subsequent runs avoid re-discovering information in
        % the subspace Y, which we already know is nearly invariant.
        % (all previous y's will be orthogonalized against on this run)
        
        clear result0 tol Iv
        % the way we re-start Arnoldi depends on whether s0 changed
        if S0(r) ~= S0(r-1)
            s0 = S0(r);
            A  = makeA(C,G,b,s0);
            result0.Vh_defl = vnext;
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
        %  as previous vectors on this run. H is pre-loaded with Up.
        result0.rho = [];
        result0.V = Y;
        n0 = size(Y,2);     % iterations will start with (n0+1)-th vector
        result0.n0 = n0;
        result0.H = Up;
        result0.exh_flag = 0;
        
        result = band_Arnoldi(n0+n(r),R0,A,tol,n0,result0);
    end
    
    if result.exh_flag > 0
        error('Cannot continue! abort on run %d',r);
    end
    
    rho = result.rho;
    H = result.H;
    if ~isempty(Up)
        H(n0+1,n0) = Up(n0+1,n0);
    end
    ph = result.Iv.ph;
    vnext = result.Vh_defl(:,ph);
    Vtmp = result.V;
    V = Vtmp(:,n0+1:result.n0);
    %  V is the set of vectors obtained from previous run of Arnoldi, not
    %  including the nearly-invariant set Y.
    
    if r < rmax
        % Extract a spanning set Y of (nearly) invariant subspace of
        %  V and possibly obtain a better expansion point based on
        %  approximate model so far
        [Y Up next_s0 v_next] = deflate(V,Y,H,vnext,s0,S0(r+1),convergence_tol);
        vnext = v_next;
        fprintf('Arnoldi run %d: %d iterations,  dim(Y)=%d\n',r,size(V,2),size(Y,2));
        S0(r+1) = next_s0;
    end
    
    Vtot = [Vtot V];   % combine them with the vectors in storage

end % loop through runs of Arnoldi

Vout = Vtot;
Yout = Y;

%% -------------------------------------------------------------

    function [Y Up_next s0_next vstart] = deflate(V,Yold,H,vnext,s0,s0_next,convergence_tol)
        ell = size(Yold,2);
        YV = [Yold V];
        ph = size(vnext,2);
        
        % determine poles, (Arnoldi relation) rel-residual error via eigs of H
        [W D] = eig(H);
        lambda = ordeig(D);
        eta = norm(vnext);
        rr = eta * abs(W(end,:)).' ./ abs(lambda); % relative residual 2-norm
		
        mu = s0 + 1./lambda;
        
        if r>1 
            rho = YV'*R0;
        end
        wt = poleweight(mu,W,YV'*c,rho,s0);
        lwt = log10(wt);
        lwt = lwt + abs(min(lwt));
        rwt = lwt/sum(lwt);  % relative weight (and log scaled) 
        
        cvrgd = rr <= convergence_tol;
        n_cvrgd =  nnz(cvrgd);
        if n_cvrgd >= min_keep
            keepers = cvrgd;
        else
            [~, rr_idx] = sort(rr);
            keepers = rr_idx(1:min_keep);
        end

        
        if isnan(s0_next)
             s0_next = autoselect_s0(mu,rr,wt,cvrgd);
            
			 % plot the poles from this run
			 figure('name',sprintf('poles r=%d, s0=%g',r,s0));
			 plot_poles(mu,rwt,rr,s0,convergence_tol);
			 title(sprintf('poles r=%d, s0=(%g,%g)',r,real(s0),imag(s0)));
			 
			 hold on
			 % and the new s0
			 plot(lscale(real(s0_next)),imag(s0_next),'gx');
			 hold off
			 drawnow
        end
        
        %  Adjust quantities for change in s0.
        if chord(s0,s0_next) < eps
            % We'd rather keep the same s0 than use a new one if the change is
            % very small.
            s0_next = s0;
            lambda_next = lambda;  % shifted lambdas
            A_next = A;
        else
            lambda_next = 1./(mu-s0_next);  % shifted lambdas
            A_next = makeA(C,G,[],s0_next);
        end
        
        % construct an orthonormal basis for nearly-invariant space to 
        %  work with the new operator (implied by the new s0) 
        %  by forming a partial Schur decomp
        [Q R] = qr(W(:,keepers),0);
        U_next = R*diag(lambda_next(keepers))/R;
        Y = YV*Q;
        
%         resid_exp = A_next(Y)-Y*U_next;
        
        Tvnext = (s0_next-s0)*A_next(vnext) + vnext;
        nrmTv = norm(Tvnext);
        vstart = Tvnext / nrmTv;
        uT = eta*nrmTv*Q(end,:);
        Up_next = [U_next; uT];
        resid_imp = vstart*uT;
        
%         chk = [mu wt rr];
        %         plot the poles from this run
        figtit = sprintf('poles r=%d, s0=%g + i%g',r,real(s0),imag(s0));
		figure('name',figtit);         
        plot_poles(mu,wt,rr,s0,convergence_tol);
        title(figtit);
        drawnow
        
        % Movie
%         figure('name',sprintf('eigmovie,r=%d',r),'Position',get(0,'Screensize'))
%         mov = eigmovie(H,vnext,s0,convergence_tol,ell,[-10.848 11.047 -3.341e+010 9.3659e+010]);
%         mov = eigmovie(H,vnext,s0,convergence_tol,ell);
%         M(midx+1:midx+length(mov)) = mov;
%         midx = midx+length(mov);
%         save('movie','M');

        
        if ~true
            %% Analysis
            % see if Arnoldi relation holds
            ar_exp = sqrt(sum(abs(A(YV) - YV*H).^2));
            ar(length(ar_exp)) = norm(vnext);
            figure('name',sprintf('arnoldi relation r=%d',r));
            h = semilogy([ar_exp(:) ar(:)],'.');
            set(h,{'Marker'},{'+','o'}');
            title(sprintf('arnoldi relation, r=%d',r))
                  
            % Compare residuals of approx eigvectors wrt A and Anext
            Z = YV*W;
            
            % Let's see how the Ritz vectors are suited to the new operator.
            rr_exp = rr_explicit(A,Z,lambda);
            rr_nadj = rr_explicit(A_next,Z,lambda_next);
            
            rr(rr<eps) = eps;
            rrs = [rr(:) rr_exp(:) rr_nadj(:)];
            rrs = sortrows(rrs,1);
            figure('name',sprintf('eig info r=%d',r))
            h = semilogy(rrs);
            set(h,{'Marker'},{'+','o','.'}.');
            line([1 length(rr)],[convergence_tol convergence_tol],'linestyle','--');
            legend('residuals',sprintf('H_%d explicit',r),sprintf('H_%d explicit',r+1),...
                'conv tol','location','best');
            title(sprintf('relative residual errors of Ritz vectors, r=%d',r));
            
            % Compare invariance of Y wrt A and Anext (see if they close)
            resid1 = A(Y) - Y*U;
            approx_resid2 = v_next*s*u_next';
                        
            AYU1 = sqrt(sum(abs(resid1).^2));
            AYU2 = sqrt(sum(abs(resid2).^2));
            AYU2_aprx = sqrt(sum(abs(approx_resid2).^2));
            nrmptb = sqrt(sum(abs(pertrb).^2));
            chk = [nrmptb(:) AYU1(:) AYU2(:) AYU2_aprx(:)];
            figure('name','schur invariance');
            h = semilogy(chk);
            set(h,{'Marker'},{'+','o','.','s'}');
            title(sprintf('Schur decomp near-invariance, r=%d',r))
            legend(sprintf('H_%d residual',r),sprintf('H_%d explicit',r),...
                sprintf('H_%d explicit',r+1),sprintf('H_%d rank-1 approx',r+1),...
                'location','best');
            
            % Check schur decomp using Up
            v = vnext/norm(vnext);
            AYUp = sqrt(sum(abs(A(Y)- [Y v]* Up).^2));
            A2YUp = sqrt(sum(abs(A_next(Y) - [Y v_next]*Up_next).^2));
            figure('name','schur2');
            chk = [AYUp(:) A2YUp(:)];
            h = semilogy(chk);set(h,{'Marker'},{'o','+'}');
            title(sprintf('near invariance relation check, r=%d',r))
            legend(sprintf('H_%d',r),sprintf('H_%d',r+1),'location','best')
        end
        
    end
end % of main function
