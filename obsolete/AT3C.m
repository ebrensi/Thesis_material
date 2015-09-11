function V_hat =  AT3C(data,m,min_ell_ess,K,samps)
% V_hat = AT3C(data,m,min_ell_ess,K,samps)
%
% m = # iterations before restart
% min_ell_ess = smallest # of restart vectors (algorithm may use more)
% K = # restarts
% samps: iterations at which poles are to be plotted
% meth: thick restart deflation step uses 0-schur vectors
%                                     or  1-orthogonalized ritz vectors
%        as the basis for the converged invariant subspace.
%
% NOTE!: This codes assume strictly complex s0  

ritz_tol = 1e-3;
tol = eps;

[C G c b] = inputdata(data);
s = getS;

if ~exist('samps','var') || isempty(samps)
    samps(1) = inf;
end
sc = 1;

ax = [-5e9 5e7 -1.5e10 1.5*2*pi*1e10];

rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
rng = rng(1:end-1) + diff(rng)/2;
S0 = sqrt(-1)*rng;

% s0 = pi*1e10;
% % s0 = (1+sqrt(-1))*pi*1e10;
% S0 = s0(ones(1,K));

s0 = S0(1);
[A r0] = makeA(C,G,b,s0);

n = length(b);
k_tot = 1;
Y = [];
V_hat = zeros(n,2*K*m);

r = 1;
while r <= K
    V = thickstart_arnoldi_basis(Y,r0,m);
    vm = size(V,2);
    V_hat(:,k_tot:k_tot+vm-1) = V;
    k_tot = k_tot+vm;
    if r < K
        [Y new_s0] = deflate(V,min_ell_ess);
        if chord(s0,new_s0) > eps
            s0 = new_s0;
            [A r0] = makeA(C,G,b,s0);
        end
    end
    r = r + 1;
end
r = r-1;

if k_tot < 2*K*m
    V_hat(:,k_tot:end) = [];
end

% This step should not be necessary but we'll use it for now
% untill we can do better.
V_hat = orth(V_hat);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
tfunc = tf_proj(V_hat,C,G,c,b);
figure;
plot_tfunc(tfunc,tfunc_real,'approx','actual');

Ck = V_hat'*C*V_hat;
Gk = V_hat'*G*V_hat;
posreal = isposreal2(Ck,Gk)

% if ~posreal
%     mu = approx_poles_proj(V_hat,C,G,c,b,s0);
%     openvar('mu');
% end

%% -------------------------------------------------------------

    function [Y new_s0] = deflate(V,min_ell)
        [mu Z rr] = approx_poles_proj(V,C,G,c,b,s0);
        cvrgd = rr < ritz_tol;
        %         wt = tf_terms_weight(tf_terms); %  weights
        %         plot_poles(mu, wt,rr,s0,ax);

        % -- pick restart vectors based on convergence (min residuals) --
        good = cvrgd;
        n_cvrgd = nnz(cvrgd);
        if n_cvrgd < min_ell
            [srr idx] = sort(rr);
            good(idx(1:min(min_ell,length(idx)))) = true;
        end

        ell = nnz(good);
        fprintf('%d converged ritz vecs, restarting with best %d\n',n_cvrgd,ell);
        new_s0 = S0(r+1);
        if ell > 0
            inv_space = Z(:,good);
            Y = make_basis([real(inv_space) imag(inv_space)]);
        else
            Y = [];
        end
    end

%%
    function [V v_next] = thickstart_arnoldi_basis(Y,v_start,m)
        V = zeros(n,m*2);
        l = size(Y,2);
        v_next = v_start/norm(v_start);

        q = real(v_start);
        V(:,1) = q/norm(q);
        kv = 1;
        if ~isreal(v_start)
            q = imag(v_start);
            q = q - (V(:,1)'*q)*V(:,1);
            V(:,2) = q/norm(q);
            kv = 2;
        end

        k = 1;
        while  k < m
            q = A(v_next);
            
            % Orthogonalize Re(q) against previous v's
            q1 = real(q);
            for i = 1:kv
                vi = V(:,i);
                q1 = q1 - (vi'*q1)*vi;
            end
            nrm_q1 = norm(q1);
            if nrm_q1 > tol
                kv = kv+1;
                V(:,kv) = q1 / nrm_q1;
            end

            if isreal(q)
                v_next = V(:,kv);
            else
                % Orthogonalize Im(q) against previous v's
                q2 = imag(q);
                for i = 1:kv
                    vi = V(:,i);
                    q2 = q2 - (vi'*q2)*vi;
                end
                nrm_q2 = norm(q2);
                if nrm_q2 > tol
                    kv = kv+1;
                    V(:,kv) = q2 / nrm_q2;
                end

                v_next = q1+sqrt(-1)*q2;
                v_next = v_next / norm(v_next);
            end
            
            % **** visualization of approximate poles ****
            if k_tot == samps(sc)
                [mu Z rr tf_terms] = approx_poles_proj(V_hat(:,1:k_tot),C,G,c,b,s0);
                [mu   rr tf_terms] = combine_large(mu, rr, tf_terms,1);

                wt = tf_terms_weight(tf_terms);
                plot_poles(mu,wt,rr,s0,ax);
                title(sprintf('r = %d, k = %d,',r,k-l));
                drawnow;
                if sc < length(samps)
                    sc = sc+1;
                end
            end
            % *****
            k = k + 1;
        end
        if kv < m*2
            V(:,kv+1:end) = [];
        end
    end
%%
%     function plottf()
%         
%     end
end
