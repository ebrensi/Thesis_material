function [Vout Hout S0 Vout2] = ROM_arnoldi(data,nvals,s0vals,convergence_tol)
% [V H] = ROM_arnoldi(datafilename,nvals,s0vals)
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
        n0 = 0;
        [A R0] = makeA(C,G,b,s0);
        result = band_Arnoldi(n(1),R0,A,tol);
        Up = [];
        Y = [];
    else
        % Subsequent runs avoid re-discovering information in
        % the subspace Y, which we already know is invariant.
        % (all previous y's will be orthogonalized against on this run)
        
        clear result0 tol Iv
        % the way we re-start Arnoldi depends on whether s0 changed
        if ~isempty(newA)
            A = newA;
            result0.Vh_defl = vnext;
            
%             result0.Vh_defl = R0;
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
    vnext = result.Vh_defl;
    Vtmp = result.V;
    V = Vtmp(:,n0+1:result.n0);
    %  V is the set of vectors obtained from previous run of Arnoldi, not
    %  including the invariant set Y.
    
    Vtot = [Vtot V];   % combine them with the vectors in storage
    
    if r < rmax
        % Extract a spanning set Y of (nearly) invariant subspace of
        %  V and possibly obtain a better expansion point based on
        %  approximate model so far
        [Vdef Y Up newA next_s0] = deflate(V,Y,H,vnext,s0,S0(r+1),convergence_tol);
        fprintf('Arnoldi run %d: %d iterations,  dim(Y)=%d\n',r,size(V,2),size(Y,2));
        s0 = next_s0;
        S0(r+1) = next_s0;
    else
        % If we're at the last run of Arnoldi, no need for invariant space.
        Vdef = deflate(V,Y,H,vnext,s0,[],convergence_tol);
    end
    Vd = [Vd Vdef];
end % loop through runs of Arnoldi


Vout = Vtot;
Vout2 = Vd;
Hout = H;  % for now just output the last H computed

%% -------------------------------------------------------------

function [Vdef Y Up newA s0_next] = deflate(V,Yold,H,vnext,s0,s0_next,convergence_tol)
   Vdef = [];  % implement this later, if at all
   
    % determine poles, (Arnoldi relation) rel-residual error via eigs of H
    [W D] = eig(H);
    lambda = diag(D);
    rr = abs(W(end,:)*norm(vnext))'./abs(lambda); % relative residual 2-norm
    mu = s0 + 1./lambda;
    
    cvrgd = rr <= convergence_tol;
    n_cvrgd =  nnz(cvrgd);
    
    if nargout > 1
        %  Autoselect a new s0 if one isn't defined
        if isinf(s0_next) || isempty(s0_next) || isnan(s0_next)
            s0_next = autoselect_s0(mu,rr,convergence_tol);
        end
        
        if chord(s0,s0_next) > eps  
            % construct new operator A if s0 is significantly different
            newA = makeA(C,G,b,s0_next);
        else
            s0_next = s0;
            newA = [];
        end
        
        if n_cvrgd > 0
            % compute Schur decomp from eig decomp (of converged eigs)
            [Q R] = qr(W(:,cvrgd),0);
            U = R*diag(lambda(cvrgd))/R;
            pertrb = norm(vnext)*Q(end,:);
            Up = [U; pertrb];
            Y = [V Yold]*Q;
        else
            Up = [];
            Y = [];
        end
    end
    
    
    
end


%%
    function s0opt = autoselect_s0(mu,rr,conv_tol)
        % Auto selects an expansion point s0 based on given convergence
        % information of poles.
        
        % eliminate poles with positive real part, negative imag part
        insignificant = real(mu)>0 | imag(mu)<0;
        mu(insignificant) = []; rr(insignificant) = [];
        
        % eliminate poles that deviate from 1 order of magnitude from most
        % poles. these are 'large' or 'small' meaning they are essentially
        % zero or Inf.
        [small large] = magvals(mu,0.7);
        outliers = small | large;
        mu(outliers) = []; rr(outliers) = [];
        
        
        % eliminate any converged poles. 
        cvrgd = rr <= conv_tol;

        
        mu(cvrgd) = []; rr(cvrgd) = [];
        
        
        
        wt = 1./abs(log10(rr));  % here weight is just inverse of rr
        wt = wt/sum(wt);
        chk1 = sortrows([mu rr wt],2)
        
        
        mwt = length(wt);
        
%         s0opt = sum(wt(1:mwt).*mu(1:mwt)); % Weighted average of poles
        
        %%Alternatively, take weighted average of logs of poles
        lr = sum(wt(1:mwt).*log10(abs(real(mu(1:mwt)))));
        li = sum(wt(1:mwt).*log10(imag(mu(1:mwt))));
        s0opt = complex(-10^lr, 10^li);
        
%         s0opt = 1i*imag(s0opt);
    end

end % of main function
