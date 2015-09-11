function [frq tfunc] = AT(data,m,min_ell_ess,K,samps)
% [frq tfunc] = AT(data,m,min_ell_ess,K,samps)

% This uses thick restart with variable l (number of restart vectors).
% This is the last version using static s0, and "hessenberg" type 
% matrix.

[Cdata G c b] = inputdata(data);
ritz_tol = 1e-5;
arnoldi_tol = eps;

if isempty(samps)
    samps(1) = inf;
end

s0 = (1+sqrt(-1))*pi*1e10;
% s0 = pi*1e10;
[A r0] = makeA(Cdata,G,b,s0);
v_start = r0/norm(r0);

[s frq]= getS();
n = length(b);

sc = 1;
k_tot = 0;
Y = []; U = []; u = 1; h_next = 1;
CCT = zeros(K*m);
CHCT = zeros(K*m);
C = eye(m);
V_hat = zeros(n,K*m);

r = 1;
while r <= K  
    [V v_next H h_next] = thickstart_arnoldi(Y,U,h_next,u,v_start,m);
    V_hat(:,(r-1)*m+1:r*m) = V(:,size(Y,2)+1:end);
    if h_next <= arnoldi_tol
        break;
    end
    
    if r < K
        [Y U X u ell] = deflate(V,H,h_next,min_ell_ess);
    end
 
    rm = r*m;
    CCT(1:rm,1:rm) = CCT(1:rm,1:rm) + C*C';
    CHCT(1:rm,1:rm) = CHCT(1:rm,1:rm) + C*H*C';
    if r < K
        CHCT(rm+1,rm) = h_next;
        C_next = [C*X zeros(r*m,m); zeros(m,ell) eye(m)];
        C = C_next;
    end
        
    v_start = v_next;
    r = r + 1;
end

if h_next <= arnoldi_tol
   fprintf('The Krylov space is maxed out.\n\n'); 
end
r = r-1;

tfunc_real = tfunc_urm(data,Cdata,G,c,b);
tfunc = transfer_function(V_hat,CCT,CHCT)';
figure;
plot_tfunc(frq,tfunc,tfunc_real,'approx','actual');

%% -------------------------------------------------------------
    function [Y U X u ell] = deflate(V,H,h_next,min_ell)
        ellm = size(V,2);
        [Z D] = eig(H);
        [S R] = schur(H);
        % In general, eigenvalues returned by eig and schur are not in 
        %  the same order, so we need to sort them.
        [ritz_vals Spos] = sort(ordeig(R));
        [ritz_vals Dpos] = sort(diag(D));
        Z = Z(:,Dpos);
        
        f = c'*V*Z;
        g = Z \ (V'*r0);
        fg = f' .* g;
        
        tf_terms = bsxfun(@rdivide, fg, 1-ritz_vals*(s-s0));
        rr = abs(h_next*Z(end,:))'./abs(ritz_vals); % relative residuals
        
        % If A is real then complex eigenvalues come in conjugate pairs,
        %  but they don't otherwise.
        Areal = isreal(s0);
        if Areal
            % combine tf_terms for conjugate pairs and hold on to only one half
            ess = imag(ritz_vals) >= 0;
            ritz_ess = ritz_vals(ess);
            cplx = imag(ritz_ess) > 0;
            tf_terms_ess = tf_terms(ess,:);
            tf_terms_ess(cplx,:) = tf_terms_ess(cplx,:) + tf_terms(~ess,:);
            tf_terms = tf_terms_ess;
            rr = rr(ess);
        end
        
        cvrgd = rr < ritz_tol;
        wt = tf_terms_weight(tf_terms); %  weights
        
        % -- pick restart vectors based on convergence (lowest residuals)
        good = cvrgd;
        n_cvrgd = nnz(cvrgd);
        if n_cvrgd < min_ell
            [srr idx] = sort(rr);
            good(idx(1:min(min_ell,length(idx)))) = true;
        end
        
        if Areal
            good_ess = good;
            good = false(ellm,1);
            good(ess) = good_ess;
            good(~ess) = good_ess(cplx);
        end
        ell = nnz(good);
        fprintf('%d converged ritz vecs, restarting with best %d\n',n_cvrgd,ell);

        goodSpos = false(ellm,1);
        goodSpos(Spos(good)) = true;
        [S R] = ordschur(S,R,goodSpos);
        
        X = S(:,1:ell);
        U = R(1:ell,1:ell);
        Y = V * X;
        u = X(end,:)';
    end
%%
    function tfunc = transfer_function(V,CCT,CHCT)
        km = length(CCT);
        bk = V'*r0;
        ck = V'*c;
        J = makeJ(V,K,m,K*m);
        rank_J = rank(J);
        if rank_J < km
            warning('MATLAB:hey','J := V''*V is rank deficient (rank %d)',rank_J);
        else
           % opts.POSDEF = true;
        end

        ckCCT = ck'*CCT;
        Fvec = V' * h_next*v_next;
        Fvec(m*(K-1)+1:K*m) = 0;
        opts.SYM = true;
        X = linsolve(J,[Fvec bk],opts);
        J_Fvec  = X(:,1);
        J_bk = X(:,2);
        B = CHCT;
        B(:,km) = B(:,km) + J_Fvec;
        
        s_s0 = s - s0;
        tfunc = zeros(length(s),1);
        for j = 1:length(s)
            tfunc(j)  = abs(ckCCT * ((CCT - s_s0(j)*B) \ J_bk));
        end
    end

%%
    function plot_poles(mu,wt,s0)
        fin = isfinite(mu);
        if ~isempty(wt)
            colormap('cool');
            scatter(real(mu(fin)),imag(mu(fin)),5,wt(fin),'filled');
            colorbar;
        else
           scatter(real(mu(fin)),imag(mu(fin)),'k.');
        end
        hold on
        plot(real(s),imag(s),'b');
        plot(real(s0),imag(s0),'r+');
        hold off
        axis(1.0e+011 *[ -8    0.5   -0.01    0.75]);
    end
%%
    function [V v_next H h_next] = thickstart_arnoldi(Y,U,h,u,v_start,m)
        l = size(Y,2);
        k = l+1;
        h_next = norm(v_start);
        ellm = m+l;
        V = zeros(n,ellm);
        H = zeros(ellm,ellm);
        if ~isempty(Y)
            V(:,1:l) = Y;
            H(1:l,1:l) = U; 
            H(l+1,1:l) = h*u';
        end
        V(:,k) = v_start;
        while  k <= ellm && h_next > arnoldi_tol
            q = A(V(:,k));
            hk = zeros(k,1);
            % Orthogonalize against previous v's
            for i = 1:k
                vi = V(:,i);
                hk(i) = vi' * q;
                q = q - hk(i)*vi;
            end
            H(1:k,k) = hk;
            h_next = norm(q);
            v_next = q / h_next;

            if k < ellm
                H(k+1,k) = h_next;
                V(:,k+1) = v_next;
            end
            
            % **** visualization of approximate poles
            k_tot = k_tot + 1;
            if k_tot == samps(sc)
                [mu rr tf_terms] = approx_poles(V(:,1:k),H(1:k,1:k),h_next);
                good = rr < ritz_tol;
                plot_poles(mu(good),wt,s,s0);
                title(sprintf('r = %d, k = %d,',r,k-l));
                drawnow;
                if sc < length(samps)
                    sc = sc+1;
                end
            end
            % *****
            k = k + 1;
        end
    end

%%
    function Jr = makeJ(V,r,k0,k_tot)
        Jr = eye(k_tot);
        if r > 1
            for j = 1:r-1
                block_j = (j-1)*k0+1:j*k0;
                after_j = j*k0+2:k_tot;
                Jr(block_j,after_j) = V(:,block_j)' * V(:,after_j);
            end
            Jr = Jr + triu(Jr,1).';
        end
    end
end
