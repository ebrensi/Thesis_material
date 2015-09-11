function [err V_hat Y_hat] =  AT4(data,m,K,ritz_tol,s0pts)
% [V_hat Y_hat] =  AT4(data,m,K,ritz_tol)
%
% m = # iterations before restart
% K = # restarts
% ritz_tol = relative residual for determining convergence in deflation
% v_type ='complex';
v_type = 'real';
if ~exist('ritz_tol','var') || isempty(ritz_tol)
    ritz_tol = 0;
end

arnoldi_tol = eps;

[C G c b] = inputdata(data);
s = getS;
tfunc_real = abs(tfunc_urm(data,C,G,c,b));

ax = [-5e9 5e7 -1.5e10 3*pi*1e10];

if exist('s0pts','var') && ~isempty(s0pts)
    S0 = pi*sqrt(-1)*s0pts*1e10;
else
    rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
    rng = rng(1:end-1) + diff(rng)/2;
    S0 = sqrt(-1)*rng;

    % s0 = pi*1e10;
    % % s0 = (1+sqrt(-1))*pi*1e10;
    % S0 = s0(ones(1,K));
end
s0 = S0(1);
[A r0] = makeA(C,G,b,s0);
v_start = r0;

n = length(b);

Y = [];
better_s0 = [];
V_hat = zeros(n,2*K*m);
Y_hat = zeros(n,2*K*m);

yk = 0;
k_tot = 0;
r = 0;

while r < K
    % generate m Arnoldi vectors, avoiding re-discovering information in
    % the subspace Y, which we already know is invariant.
    if exist('v_type','var') && strcmp(v_type,'real')
        [V v_next] = thickReal4(A,Y,v_start,m);
    else
        v_type = 'complex';
        [V v_next] = thickCplx(A,Y,v_start,m);
        %         [V v_next] = arnoldiT(A,Y,v_start,m);
    end
    vm = size(V,2);

    % combine them with Arnoldi vectors in storage
    V_hat(:,k_tot+1:k_tot+vm) = V;
    k_tot = k_tot+vm;

    r = r+1;
    if r < K
        % Extract a (nearly) invariant subspace Y from V
        if ritz_tol
            [Y better_s0] = deflate(V,false);
        end
        l = size(Y,2);

        if l > 0
            % add Y to the invariant subspace already discovered
            Y_hat(:,yk+1:yk+l) = Y;
            yk = yk+l;
        end
        % all previous y's will be orthogonalized against on the next run
        Y = Y_hat(:,1:yk);

        % If deflate didn't return a better s0 then just use the next one
        % on the list.
        if ~isempty(better_s0)
            next_s0 = better_s0;
        else
            next_s0 = S0(r+1);
        end

        if chord(s0,next_s0) > eps
            s0 = next_s0;
            [A r0] = makeA(C,G,b,s0);
            v_start = r0;
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


if strcmp(v_type,'complex')
    V_hat = orth([real(V_hat) imag(V_hat)]);
end

tfunc = tf_proj(V_hat,C,G,c,b);
err = plot_tfunc(tfunc,tfunc_real,'approx','actual');
title(sprintf('%s (%s):  m=%d, r=%d, err=%g','AT4',v_type,m,K,err));

Ck = V_hat'*C*V_hat;
Gk = V_hat'*G*V_hat;
posreal = isposreal2(Ck,Gk)

% if ~posreal
%     mu = approx_poles_proj(V_hat,C,G,c,b,s0);
%     format long g
%     mu(real(mu)>0)
% %     openvar('mu');
% end

%% -------------------------------------------------------------

    function [Y better_s0] = deflate(V,q)
        [mu Z rr tf_terms] = approx_poles_proj(V,C,G,c,b,s0);
        cvrgd = rr < ritz_tol;


        %         wt = tf_terms_weight(tf_terms); %  weights
        %         wt = wt./sum(wt);
        %         %         chk_wt = sortrows([mu wt rr cvrgd],2);
        %         %         chk_rr = sortrows(chk_wt,3);
        %         plot_poles(mu, wt,rr,s0,ax);
        %         figure;
        % %         tfunc = tf_proj(V_hat(:,1:k_tot),C,G,c,b);
        %         tfunc = tf_proj(V,C,G,c,b);
        %         plot_tfunc(tfunc,tfunc_real,'approx','actual');
        %         figure;
        %
        % -- pick restart vectors based on convergence (min residuals)
        good = cvrgd;
        n_cvrgd = nnz(cvrgd);

        ell = nnz(good);
        if ~exist('q','var') || ~q
            fprintf('%d converged ritz vecs, restarting with best %d\n',n_cvrgd,ell);
        end
        better_s0 = [];

        if ell > 0
            inv_space = Z(:,good);
            if isreal(V);
                Y = orth([real(inv_space) imag(inv_space)]);
            else
                Y = orth(inv_space);
            end
        else
            Y = [];
        end
    end
%%
    function [V v_next] = thickCplx(A,Y,v_start,m)
        V = zeros(n,m);
        v_start = orthagainst(v_start,Y);
        V(:,1) = v_start/norm(v_start);
        k = 1;
        while  k <= m
            q = A(V(:,k));
            q = orthagainst(q,[Y V(:,1:k)]);
            v_next = q / norm(q);

            k = k + 1;
            if k <= m
                V(:,k) = v_next;
            end
        end
    end

%%
    function [V v_next] = thickReal(A,Y,v_start,m)
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
    function [V v_next] = thickReal2(A,Y,v_start,m)
        V = zeros(n,m*2);
        k = 0;
        kv = 0;
        while  k < m
            if k
                q = A(v_next);
            else
                q = v_start;
            end
            % Orthogonalize [Re(q) Im(q)] against Y and previous v's
            q = orthagainst([real(q) imag(q)],[Y V(:,1:kv)]);

            q1 = q(:,1);
            if norm(q1) > arnoldi_tol
                kv = kv+1;
                V(:,kv) = q1;
            end

            q2 = q(:,2);
            if norm(q2) > arnoldi_tol
                kv = kv+1;
                V(:,kv) = q2;
            end

            v_next = q1+sqrt(-1)*q2;
            v_next = v_next / norm(v_next);

            k = k + 1;
        end

        if kv < m*2
            V(:,kv+1:end) = [];
        end
    end
%%
% this version only stores real part of Av
    function [V v_next] = thickReal3(A,Y,v_start,m)
        V = zeros(n,m);
        k = 0;
        exhausted = false;
        while  k < m && ~exhausted
            if k
                q = A(v_next);
            else
                q = v_start;
            end
            % Orthogonalize Re(q) against Y and previous v's
            q1 = orthagainst(real(q),[Y V(:,1:k)]);

            nrm_q1 = norm(q1);
            if nrm_q1 > arnoldi_tol
                v_next = q1 / nrm_q1;
            else
                % if Re(q) is linearly dependent, try Im(q)
                q2 = orthagainst(imag(q),[Y V(:,1:k) ]);
                nrm_q2 = norm(q2);
                if nrm_q2 > arnoldi_tol
                    v_next = q2 / nrm_q2;
                else
                    exhausted = true;
                end
            end

            if ~exhausted
                k = k + 1;
                V(:,k) = v_next;
            else
                fprintf('exhausted at k = %d',k);
            end
        end
        if k < m
            V(:,k+1:end) = [];
        end
    end

%%
   % this version computes and stores Re and Im parts, but only passes on
   % Real part.
    function [V v_next] = thickReal4(A,Y,v_start,m)
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
                
                v_next = V(:,kv);
            end
            k = k + 1;
        end
        if kv < m*2
            V(:,kv+1:end) = [];
        end
    end

%%
end % of main function
