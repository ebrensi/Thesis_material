function [err V_hat Y_hat] =  AT4PG(data,m,min_ell,K,meth)
% [V_hat Y_hat] =  AT4(data,m,min_ell,K,meth)
%
% m = # iterations before restart
% min_ell = smallest # of restart vectors (algorithm may use more)
% K = # restarts
% samps: iterations at which poles are to be plotted
% meth: thick restart deflation step uses 0-schur vectors
%                                     or  1-orthogonalized ritz vectors
%        as the basis for the converged invariant subspace.

%   AT4PG is an experiment with deflating out "pretty good" vectors,
%   which are not quite good enough to be considered invariant but whose
%   elimination from the search may be desirable. These vectors will be 
%   orthogonalized against on the next run of Arnoldi but not kept around.  
% (and it doesn't really improve anything)

ritz_tol = 1e-4;  % acceptable relative error to consider a pole converged
arnoldi_tol = eps; 

[C G c b] = inputdata(data);
s = getS;

ax = [-5e9 5e7 -1.5e10 3*pi*1e10];

rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
rng = rng(1:end-1) + diff(rng)/2;
S0 = sqrt(-1)*rng;

% s0 = pi*1e10;
% % s0 = (1+sqrt(-1))*pi*1e10;
% S0 = s0(ones(1,K));

s0 = S0(1);
[A r0] = makeA(C,G,b,s0);
v_start = r0;

n = length(b);

Y = [];
almostY = [];
V_hat = zeros(n,2*K*m);
Y_hat = zeros(n,2*K*m);

yk = 0;
k_tot = 0;
r = 0;

while r < K
    % generate m Arnoldi vectors, avoiding re-discovering information in
    % the subspace Y, which we already know is invariant.
    if strcmp(meth,'real')
        [V v_next] = thickReal(Y,almostY,v_start,m);
    else
        meth = 'complex';
        [V v_next] = thickCplx(Y,almostY,v_start,m);
    end
    vm = size(V,2);

    % combine them with Arnoldi vectors in storage
    V_hat(:,k_tot+1:k_tot+vm) = V;
    k_tot = k_tot+vm;

    r = r+1;
    if r < K
        % Extract a (nearly) invariant subspace Y from V
        %  i.e. distill V into the eigen-information it contains
        [Y almostY better_s0] = deflate(V,min_ell);
        l = size(Y,2);

        if l > 0
            % add Y to the invariant subspace already discovered
            Y_hat(:,yk+1:yk+l) = Y;
            yk = yk+l;
        end
        % all previous y's will be orthogonalized against on the next run
        Y = Y_hat(:,1:yk);
        
        if isempty(better_s0)
            better_s0 = S0(r+1);
        end
        if chord(s0,better_s0) > eps
            s0 = better_s0;
            [A r0] = makeA(C,G,b,s0);
            v_start = r0;  % this may not be good
%             v_start = v_next;
        else
            v_start = v_next;
        end
    end
end

if k_tot < 2*K*m
    V_hat(:,k_tot+1:end) = [];
end
if yk < 2*K*m
    Y_hat(:,yk+1:end) = [];
end

if ~isreal(V_hat)
%     V_hat = orth([real(V_hat) imag(V_hat)]);
    V_hat = make_basis([real(V_hat) imag(V_hat)]);
end

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
tfunc = tf_proj(V_hat,C,G,c,b);
err = plot_tfunc(tfunc,tfunc_real,'approx','actual');
title(sprintf('%s (%s):  m=%d, l= %d, r=%d, err=%g','AT4',meth,m,min_ell,K,err));
Ck = V_hat'*C*V_hat;
Gk = V_hat'*G*V_hat;
posreal = isposreal2(Ck,Gk)

if ~posreal
    mu = approx_poles_proj(V_hat,C,G,c,b,s0);
    format long g 
    mu(real(mu)>0)
%     openvar('mu');
end

%% -------------------------------------------------------------

    function [Y almostY better_s0] = deflate(V,min_ell)
        [mu Z rr tf_terms] = approx_poles_proj(V,C,G,c,b,s0);
        cvrgd = rr < ritz_tol;
        wt = tf_terms_weight(tf_terms); %  weights
        plot_poles(mu, wt,rr,s0,ax);

        % -- pick restart vectors based on convergence (min residuals)
        n_cvrgd = nnz(cvrgd);
        pretty_good = false(size(mu));
        
        % If the user wants to restart with more vectors than have 
        % converged, we put some "pretty good" ones in too, but 
        % those will be treated differently. (not permanent)
        if n_cvrgd < min_ell
            [srr idx] = sort(rr);
            pretty_good(idx(n_cvrgd+1:min(min_ell,length(idx)))) = true;
        end

        n_ok = nnz(pretty_good);
        better_s0 = [];

        if n_cvrgd > 0
            inv_space = Z(:,cvrgd);
            if isreal(V)
                Y = orth([real(inv_space) imag(inv_space)]);
            else
                Y = orth(inv_space);                
            end
        else
            Y = [];
        end
        
        if n_ok > 0
            ok_space = Z(:,pretty_good);
            if isreal(V)
                almostY = orth([real(ok_space) imag(ok_space)]);
            else
                almostY = orth(ok_space);                
            end
        else
            almostY = [];
        end
        
        fprintf('%d converged, restarting with best %d\n',n_cvrgd,n_cvrgd+n_ok);
    end
%%
    function [V v_next] = thickCplx(Y,almostY,v_start,m)
        l = size(Y,2);
        V = zeros(n,m);
%         v_start2 = orthagainst(v_start,Y);
%         norm(v_start-v_start2)/norm(v_start)
        V(:,1) = v_start/norm(v_start);
        k = 1;
        while  k <= m
            q = A(V(:,k));
            q = orthagainst(q,[Y almostY V(:,1:k)]);
            v_next = q / norm(q);

            k = k + 1;
            if k <= m
                V(:,k) = v_next;
            end
        end
    end

%%
    function [V v_next] = thickReal(Y,v_start,m)
        V = zeros(n,m*2);
        k = 0;
        kv = 0;
        while  k < m
            if k
                q = A(v_next);
            else
                q = v_start;
            end
            % Orthogonalize Re(q) against Y and previous v's
            q1 = orthagainst(real(q),[Y V(:,1:kv)]);

            nrm_q1 = norm(q1);
            if nrm_q1 > arnoldi_tol
                kv = kv+1;
                V(:,kv) = q1 / nrm_q1;
            end

            if isreal(q)
                v_next = V(:,kv);
            else
                % Orthogonalize Im(q) against Y and previous v's
                q2 = orthagainst(imag(q),[Y V(:,1:kv) ]);
                nrm_q2 = norm(q2);
                if nrm_q2 > arnoldi_tol
                    kv = kv+1;
                    V(:,kv) = q2 / nrm_q2;
                end

                v_next = q1+sqrt(-1)*q2;
                v_next = v_next / norm(v_next);
            end
            k = k + 1;
        end
        if kv < m*2
            V(:,kv+1:end) = [];
        end
    end
%%
    function [V v_next] = thickReal2(Y,v_start,m)
        V = zeros(n,m*2);
        k = 0;
        kv = 0;
        while  k < m
            if k
                q = A(v_next);
            else
                q = v_start;
            end
            % Orthogonalize Re(q) and Im(q) against Y and previous v's
            q12 = orthagainst([real(q) imag(q)],[Y V(:,1:kv)]);
            
            q1 = q12(1);
            nrm_q1 = norm(q1);
            if nrm_q1 > arnoldi_tol
                kv = kv+1;
                V(:,kv) = q1 / nrm_q1;
            end

            if isreal(q)
                v_next = V(:,kv);
            else
                % Orthogonalize Im(q) against previous v's
                q2 = orthagainst(q12(2),V(:,kv));
                nrm_q2 = norm(q2);
                if nrm_q2 > arnoldi_tol
                    kv = kv+1;
                    V(:,kv) = q2 / nrm_q2;
                end

                v_next = q1+sqrt(-1)*q2;
                v_next = v_next / norm(v_next);
            end
            k = k + 1;
        end
        if kv < m*2
            V(:,kv+1:end) = [];
        end
    end
%%
    function Xout = orthagainst(Xin,basis)
        for i = 1:size(basis,2)
            bi = basis(:,i);
            Xin = Xin - bsxfun(@times,(bi'*Xin),bi);
        end
        Xout = Xin;
    end
end
