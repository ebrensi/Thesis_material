function V_hat =  AT3(data,m,min_ell_ess,K)
% [frq tfunc] = AT3(data,m,min_ell_ess,K)
%
% m = # iterations before restart
% min_ell_ess = smallest # of restart vectors (algorithm may use more)
% K = # restarts
% samps: iterations at which poles are to be plotted
% meth: thick restart deflation step uses 0-schur vectors
%                                     or  1-orthogonalized ritz vectors
%        as the basis for the converged invariant subspace.

ritz_tol = 1e-3;

[C G c b] = inputdata(data);
s = getS;

ax = [-5e9 5e7 -1.5e10 1.5*2*pi*1e10];

rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
rng = rng(1:end-1) + diff(rng)/2;
S0 = sqrt(-1)*rng;

% s0 = pi*1e10;
% % s0 = (1+sqrt(-1))*pi*1e10;
% S0 = s0(ones(1,K));

s0 = S0(1);
[A r0] = makeA(C,G,b,s0);
v_start = r0/norm(r0);

n = length(b);
k_tot = 0;
Y = [];
V_hat = zeros(n,K*m);

r = 1;
while r <= K
    [V v_next] = thickstart_arnoldi_basis(Y,v_start,m);
    V_hat(:,(r-1)*m+1:r*m) = V;
    
    if r < K
        [Y new_s0] = deflate(V,min_ell_ess);
        if chord(s0,new_s0) > eps
            s0 = new_s0;
            [A r0] = makeA(C,G,b,s0);
            v_start = r0;
        end
    else
        v_start = v_next;
    end
    
    r = r + 1;
end

r = r-1;

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
tfunc = tf_proj(V_hat,C,G,c,b);
figure;
plot_tfunc(tfunc,tfunc_real,'approx','actual');

Vreal = make_basis([real(V_hat) imag(V_hat)]);
Ck = Vreal'*C*Vreal;
Gk = Vreal'*G*Vreal;
posreal = isposreal2(Ck,Gk)


%% -------------------------------------------------------------

    function [Y new_s0] = deflate(V,min_ell)
        [mu Z rr] = approx_poles_proj(V,C,G,c,b,s0);
        cvrgd = rr < ritz_tol;
        %         wt = tf_terms_weight(tf_terms); %  weights
        %         plot_poles(mu, wt,rr,s0,ax);

        % -- pick restart vectors based on convergence (min residuals)
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
            Y = orth(Z(:,good));
        else
            Y = [];
        end
    end
%%
    function [V v_next] = thickstart_arnoldi_basis(Y,v_start,m)
        l = size(Y,2);
        V = zeros(n,m);
        V(:,1) = v_start/norm(v_start);
        k = 1;
        while  k <= m
            q = A(V(:,k));
            if l > 0
                % Orthogonalize against y's
                for i = 1:l
                    yi = Y(:,i);
                    q = q - (yi'*q)*yi;
                end
            end

            % Orthogonalize against previous v's
            for i = 1:k
                vi = V(:,i);
                q = q - (vi'*q)*vi;
            end
            v_next = q / norm(q);

            if k < m
                V(:,k+1) = v_next;
            end

            k = k + 1;
        end
    end
%%
end
