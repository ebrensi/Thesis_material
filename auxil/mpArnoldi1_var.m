function [Vout Yout S0] = mpArnoldi1_var(data,nvals,s0vals,convergence_tol)
% [Vout Yout S0] = mpArnoldi1(data,nvals,s0vals,convergence_tol)
%
% nvals = # iterations for run(s)
% s0vals = (optional) s0 points for each run
%
% example:  ROM_arnoldi('ex308s_1',[50 10 10],(1+1.5i)*pi*1e10)
%
%  this does three runs of thick-restart band_Arnoldi of 50, 10 and 10
%  iterations, starting with a specified initial s0 point and letting
%  the adaptive routine choose the other two expansion points.
%
% example: ROM_arnoldi('ex308s_1',[50 inf 10],[1+1.5i, 1i, 3+2i]*pi*1e10)
%
%  Three runs of thick restart band_arnoldi with the first and third s0
%   predetermined and the second one to be chosen by the adaptive algorithm.

% in this version we include Y as part of a new start block R0 = [Y vnext],
%  rather than pre-loading V with it.

%% Internally defined parameters
defl_tol = sqrt(eps);  % deflation tolerance for band_Arnoldi
s0init = (0.2+5i)*2*pi*1e9; % default initial s0 (should be generally located)

%% Setup
% create vector of s0 values.  Those to be autoselected are indicated by inf.
rmax = size(nvals,2);
S0 = inf(1,rmax);
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
Vd = [];

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
        [A R0] = makeA(C,G,b,s0);
        result = band_Arnoldi(n(r),R0,A,tol);
        Y = [];
        ell = 0;
    else
        % Subsequent runs avoid re-discovering information in
        % the subspace Y, which we already know is invariant.
        % (all previous y's will be orthogonalized against on this run)
        
        % the way we re-start Arnoldi depends on whether s0 changed
        if S0(r) ~= S0(r-1)
            s0 = S0(r);
            A = makeA(C,G,b,s0);
            
            tol.defl_flag = 1;
            tol.defl_tol = defl_tol;
            tol.normA_flag = 1;
            tol.normA = 0;
        else
            tol = result.tol;
        end
        
        ell = size(Y,2);
        R0 = [Y vnext];
        result = band_Arnoldi(n(r)+ell,R0,A,tol);
    end
    
    if result.exh_flag > 0
        error('Cannot continue! abort on run %d',r);
    end
    
    H = result.H;
    V = result.V;
    ph = result.Iv.ph;
    vnext = result.Vh_defl(:,ph);
    
    if r < rmax
        % Extract a spanning set Y of (nearly) invariant subspace of
        %  V and possibly obtain a better expansion point based on
        %  approximate model so far
        [Y next_s0] = deflate(V,H,vnext,s0,S0(r+1),convergence_tol);
        H = result.H(ell+1:end,ell+1:end);
        V = result.V(:,ell+1:end);
        fprintf('Arnoldi run %d: %d iterations,  dim(Y)=%d\n',r,size(V,2),size(Y,2));
        S0(r+1) = next_s0;
    end
    
    Vtot = [Vtot V];   % combine them with the vectors in storage
end % loop through runs of Arnoldi


Vout = Vtot;
Yout = Y;

%% -----------------------------------------------------------------------------

    function [Y s0_next] = deflate(V,H,vnext,s0,s0_next,convergence_tol)
        % see if Arnoldi relation holds 
        ar_exp = sqrt(sum(abs(A(V) - V*H).^2));
        ar(length(ar_exp)) = norm(vnext);
        figure('name',sprintf('arnoldi relation r=%d',r));
        h = semilogy([ar_exp(:) ar(:)],'.');
        set(h,{'Marker'},{'+','o'}');
        title(sprintf('arnoldi relation, r=%d',r))
        
        % determine poles, (Arnoldi relation) rel-residual error via eigs of H
        [W D] = eig(H);
        lambda = ordeig(D);
        eta = norm(vnext);
        rr = eta * abs(W(end,:)).' ./ abs(lambda); % relative residual 2-norm
        rr(rr<eps) = eps;
        
        
        % sort eigs by magnitude  (we do this to match the schur decomp)
        [lambda mag_idx] = sort(lambda);
        W = W(:,mag_idx);
        rr = rr(mag_idx);

        mu = s0 + 1./lambda;
        cvrgd = rr <= convergence_tol;
        n_cvrgd =  nnz(cvrgd);
        
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
        
        if n_cvrgd
            % construct a nearly-invariant space to work with the new
            % operator (implied by the new s0)
            % via partial Schur decomp 
            
%             [Q R] = qr(W(:,cvrgd),0);
%             U = R*diag(lambda(cvrgd))/R;
%             U_next = R*diag(lambda_next(cvrgd))/R;
%             uT = Q(end,:);  % the last row of Q
%             Up = [U_next; eta*uT];
            Y = V*W(:,cvrgd);
        else
            Y = [];
        end
        
        %% Analysis
        % Compare residuals of approx eigvectors wrt A and Anext
        Z = V*W;
        
        % Let's see how the Ritz vectors are suited to the new operator.
        rr_exp = rr_explicit(A,Z,lambda);
        rr_nadj = rr_explicit(A_next,Z,lambda_next);
        
        rrs = [rr(:) rr_exp(:) rr_nadj(:)];
        rrs = sortrows(rrs,1);
        figure('name',sprintf('eig info r=%d',r))
        h = semilogy(rrs);
        set(h,{'Marker'},{'+','o','.'}.');
        line([1 length(rr)],[convergence_tol convergence_tol],'linestyle','--');
        legend('residuals',sprintf('H_%d explicit',r),sprintf('H_%d explicit',r+1),'conv tol','location','best');
        title(sprintf('relative residual errors of Ritz vectors, r=%d',r));
    end
end % of main function
