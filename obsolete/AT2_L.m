function [frq tfunc] = AT2_L(data,m,min_ell_ess,K,samps,matdef)
% [frq tfunc] = AT2(data,m,min_ell_ess,K,samps,matdef)
% experiments with elimination of linear dependence in V.

[Cdata G c b] = inputdata(data);

ritz_tol = 1e-5;
arnoldi_tol = eps;

if isempty(samps)
    samps(1) = inf;
end
sc = 1;

s = getS;

% rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
% rng = rng(1:end-1) + diff(rng)/2;
% S0 = 1e3 + sqrt(-1)*rng;
% s0 = S0(1);

s0 = pi*1e10;
% s0 = (1+sqrt(-1))*pi*1e10;
S0 = s0(ones(1,K));

[A r0] = makeA(Cdata,G,b,s0);
v_start = r0/norm(r0);

n = length(b);
k_tot = 0;
Y = []; U = []; u = 1; h_next = 1;
V_hat = zeros(n,K*m);

r = 1;
vcount = 0;
while r <= K  
    [V v_next H h_next] = thickstart_arnoldi(Y,U,h_next,u,v_start,m);
    
    if matdef && r > 1
        [V_nice bad] = matDeflate([V_hat(:,1:vcount) V(:,size(Y,2)+1:end)]);
        vcount = size(V_nice,2);
        fprintf('\nV deflate: elminated %d vectors, dim(V_hat) = %d\n',bad,vcount);
        V_hat(:,1:vcount) = V_nice;
    else
        V_hat(:,vcount+1:vcount+m) = V(:,size(Y,2)+1:end);
        vcount = vcount + m;
    end
    
    if h_next <= arnoldi_tol
        break;
    end
    
    if r < K
        [Y U X u ell new_s0] = e_deflate(V,H,h_next,min_ell_ess);
        if chord(s0,new_s0) > eps
            s0 = new_s0;
            A = makeA(Cdata,G,s0);
        end
    end

    v_start = v_next;
    r = r + 1;
end

if h_next <= arnoldi_tol
   fprintf('The Krylov space is maxed out.\n\n'); 
end
r = r-1;

tfunc_real = abs(tfunc_urm(data,Cdata,G,c,b));
tfunc = tf_proj(V_hat,Cdata,G,b,c);
figure;
plot_tfunc(frq,tfunc,tfunc_real,'approx','actual');


%% -------------------------------------------------------------
    function [Y U X u ell new_s0] = e_deflate(V,H,h_next,min_ell)
        % note: we don't really need to compute weights except for analysis
        %  purposes.  
        ellm = size(V,2);
        [Z D] = eig(H);
        ritz_vals = diag(D);
       
        rr = abs(h_next*Z(end,:))'./abs(ritz_vals); % relative residuals
        
        % If A is real then complex eigenvalues come in conjugate pairs,
        %  but they don't otherwise.
        Areal = isreal(s0);
        if Areal
            % combine tf_terms for conjugate pairs and hold on to only one half
            ess = imag(ritz_vals) >= 0;
            ritz_ess = ritz_vals(ess);
            cplx = imag(ritz_ess) > 0;
            rr = rr(ess);
%         else
%             % if A is complex, eigs don't come in complex pairs but since 
%             %  poles of the pade approx do, we get extra ritz pairs.
%             extra_ritz_vals = conj(ritz_vals) ./ (1-2*sqrt(-1)*imag(s0).*conj(ritz_vals));
%             rr2 = abs(h_next*conj(Z(end,:)))'./abs(extra_ritz_vals); % relative residuals
        end
       
        cvrgd = rr < ritz_tol;
        
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
                plot_poles(mu(good),wt,rr,s0);
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
    function [Vgood n_bad] = matDeflate(Vbad)
        [UU,SS,VV] = svd(Vbad);
        sing_vals = diag(SS);
        good = sing_vals > eps / 2;
        n_bad = nnz(~good);
        Vgood = UU(:,good);       
    end
end
