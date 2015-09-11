function [frq tfunc] = AT2(data,m,min_ell_ess,K,samps,meth)
% [frq tfunc] = AT2(data,m,min_ell_ess,K,samps)
%
% m = # iterations before restart
% min_ell_ess = smallest # of restart vectors (algorithm may use more)
% K = # restarts
% samps: iterations at which poles are to be plotted
% meth: thick restart deflation step uses 0-schur vectors 
%                                     or  1-orthogonalized ritz vectors
%        as the basis for the converged invariant subspace.  

[C G c b] = inputdata(data);

ritz_tol = 1e-5;
arnoldi_tol = eps;

if isempty(samps)
    samps(1) = inf;
end
sc = 1;

frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;

rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
rng = rng(1:end-1) + diff(rng)/2;
S0 = 1e3 + sqrt(-1)*rng;

% s0 = pi*1e10;
% % % s0 = (1+sqrt(-1))*pi*1e10;
% S0 = s0(ones(1,K));

s0 = S0(1);
[A r0] = makeA(C,G,b,s0);
v_start = r0/norm(r0);

n = length(b);
k_tot = 0;
Y = []; U = []; u = 1; h_next = 1;
V_hat = zeros(n,K*m);

r = 1;
while r <= K  
    [V v_next H h_next] = thickstart_arnoldi(Y,U,h_next,u,v_start,m);
    V_hat(:,(r-1)*m+1:r*m) = V(:,size(Y,2)+1:end);
   
    if h_next <= arnoldi_tol
        break;
    end
    
    if r < K
        if meth == 1 
            [Y U X u ell new_s0] = e_deflate(V,H,h_next,min_ell_ess);
        else
            [Y U X u ell new_s0] = s_deflate(V,H,h_next,min_ell_ess);
        end

        if chord(s0,new_s0) > eps
            s0 = new_s0;
            A = makeA(C,G,b,s0);
        end
    end

%     tfunc = transfer_function_alt(V_hat(:,1:k_tot),C,G,c,b,s);
%     figure;
%     plot_tfunc(tfunc);
%     drawnow;

    v_start = v_next;
    r = r + 1;
end

if h_next <= arnoldi_tol
   fprintf('The Krylov space is maxed out.\n\n'); 
end
r = r-1;

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
tfunc = tf_proj(V_hat,C,G,c,b,s);
figure;
plot_tfunc(frq,tfunc,tfunc_real,'approx','actual');

Ck = V_hat'*C*V_hat;
Gk = V_hat'*G*V_hat;
bk = V_hat'*b;
ck = V_hat'*c;
posreal = isposreal(Ck,Gk,ck,bk,s0)


%% -------------------------------------------------------------
    function [Y U X u ell new_s0] = e_deflate(V,H,h_next,min_ell)
        % note: we don't really need to compute weights except for analysis
        %  purposes.  
        ellm = size(V,2);
        [Z D] = eig(H);
        ritz_vals = diag(D);
        
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
%         else
%             % if A is complex, eigs don't come in complex pairs but since 
%             %  poles of the pade approx do, we get extra ritz pairs.
%             extra_ritz_vals = conj(ritz_vals) ./ (1-2*sqrt(-1)*imag(s0).*conj(ritz_vals));
%             rr2 = abs(h_next*conj(Z(end,:)))'./abs(extra_ritz_vals); % relative residuals
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
        new_s0 = S0(r+1);
        % shift the ritz values for the new s0
        ds0 = new_s0 - s0;
%         new_ritz_vals = ritz_vals(good) ./ (1-ds0*ritz_vals(good));
        
        X = orth(Z(:,good));
        U = diag(ritz_vals(good));
        U = U*inv(eye(ell)-ds0*U);
        Y = V * X;
        u = X(end,:)';
    end

%%
    function [Y U X u ell new_s0] = s_deflate(V,H,h_next,min_ell)
        % note: we don't really need to compute weights except for analysis
        %  purposes.  
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
%         else
%             % if A is complex, eigs don't come in complex pairs but since 
%             %  poles of the pade approx do, we get extra ritz pairs.
%             extra_ritz_vals = conj(ritz_vals) ./ (1-2*sqrt(-1)*imag(s0).*conj(ritz_vals));
%             rr2 = abs(h_next*conj(Z(end,:)))'./abs(extra_ritz_vals); % relative residuals
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
              
%         gg = [ritz_vals(good) rr(good) Spos(good)];
        goodSpos = false(ellm,1);
        goodSpos(Spos(good)) = true;
        [S R] = ordschur(S,R,goodSpos);
        
        X = S(:,1:ell);
        U = R(1:ell,1:ell);
                
        % Since s0 is different on next run, A(s0) will be different too
        %  so we need to transform this (nearly) A invariant space to 
        %  be invariant wrt to the new A.
        new_s0 = S0(r+1);
        ds0 = new_s0 - s0;
%         new_ritz_vals = ritz_vals(good) ./ (1-ds0*ritz_vals(good));
        I = eye(ell);
%         X = orth(X*inv(I-ds0*U));
        U = U*inv(I-ds0*U);
        
        Y = V * X;
        u = X(end,:)';
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
                [mu rr tf_terms] = approx_poles(V(:,1:k),H(1:k,1:k),h_next,c,r0,s,s0);
                plot_poles(mu,wt,rr,s0);
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
end
