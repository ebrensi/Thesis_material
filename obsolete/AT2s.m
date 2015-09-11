function [frq tfunc] = AT2s(data,m,min_ell_ess,K,samps)
% [frq tfunc] = AT2s(data,m,min_ell_ess,K,samps)

% 10/26/2008

C = []; G = []; b = []; c = [];
load(data);
Cdata = C;

ritz_tol = 1e-5;
arnoldi_tol = eps;

if isempty(samps)
    samps(1) = inf;
end


frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;

rng = linspace(0,max(imag(s))+min(imag(s)),K+1);
rng = rng(1:end-1) + diff(rng)/2;
S0 = 1 + sqrt(-1)*rng;
s0 = S0(1);

% s0 = pi*1e10;
% % s0 = (1+sqrt(-1))*pi*1e10;
% S0 = s0(ones(K,1));

C_s0G = Cdata-s0*G;
r0 = -C_s0G\b;
v_start = r0/norm(r0);

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
        [Y U X u ell new_s0] = deflate(V,H,h_next,min_ell_ess);
        if chord(s0,new_s0) > eps
            s0 = new_s0;
            C_s0G = Cdata-s0*G;
            r0 = -C_s0G\b;
            v_start = r0/norm(r0);
        end
    end
 
    rm = r*m;
    CCT(1:rm,1:rm) = CCT(1:rm,1:rm) + C*C';
    CHCT(1:rm,1:rm) = CHCT(1:rm,1:rm) + C*H*C';
    if r < K
        CHCT(rm+1,rm) = h_next;
        C_next = [C*X zeros(r*m,m); zeros(m,ell) eye(m)];
        C = C_next;
    end

%     tfunc = transfer_function_alt(V_hat(:,1:k_tot),Cdata,G,b,c,s);
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

tfunc_real = transfer_function_real(Cdata,G,b,c,s);

% tfunc = transfer_function(V_hat,CCT,CHCT);
% figure;
% plot_tfunc(tfunc);

tfunc = transfer_function_alt(V_hat,Cdata,G,b,c,s);
figure;
plot_tfunc(tfunc,tfunc_real,'approx','actual');

%% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end
%%
    function [Y U X u ell new_s0] = deflate(V,H,h_next,min_ell)
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
        
%         E = A(V*S) - V*S*R;
%         normsr = sqrt(sum(abs(R).^2))';
%         ee = sqrt(sum(abs(E).^2))' ./ normsr;
%         rr2 = abs(h_next*S(end,:))'./normsr;
%         rr2 = rr2(Spos);
        
              
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
        wt = max(abs(tf_terms),[],2); %  weights
        
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
        X = orth(X*inv(I-ds0*U));
        U = U*inv(I+ds0*U);
        
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
    function tfunc = transfer_function_alt(Vk,C,G,c,b,s)
        Ls = length(s);
        Ck = Vk'*C*Vk;
        Gk = Vk'*G*Vk;
        bk = Vk'*b;
        ck = Vk'*c;
        tfunc = zeros(1,Ls);
        for j = 1:Ls
            tfunc(j) = ck'*((s(j)*Gk-Ck)\bk);
        end
    end

%%
    function tfunc_real = transfer_function_real(C,G,c,b,s)
        exname = data(3:end);
        tfname = sprintf('tfunc_%s.mat',exname);
        if exist(tfname,'file') == 2
            vars = load(tfname);
            tfunc_real = vars.tfunc_real;
        else
            Ls = length(s);
            tfunc_real = zeros(1,Ls);
            for j = 1:length(s)
                tfunc_real(j) = c'*((s(j)*G-C)\b);
            end
            save(tfname,'tfunc_real','frq');
        end
    end

%%
    function plot_tfunc(tfunc,alt_tfunc,name1,name2)
        Tfunc = abs(tfunc);
        % plot it on a log scale
        if exist('alt_tfunc','var')
            alt_Tfunc = abs(alt_tfunc);
            loglog(frq,alt_Tfunc,'r--',frq,Tfunc);
            err = norm(alt_Tfunc-Tfunc,inf) / norm(alt_Tfunc,inf);
            legend(name2, sprintf('%s (rel diff=%.4g)',name1,err) );
        else
            loglog(frq,Tfunc);
        end
        xlabel('frq','fontsize',12);
        ylabel('|H_K(s)|','fontsize',12,'Rotation',90);
        title(sprintf('m = %d, K = %d, l_{min} = %d',m,K,min_ell_ess));
        drawnow
        %         saveas(h,sprintf('A1_%s_%d.png',data,k));
        %         close h
    end


%%
    function plot_poles(mu,wt,s,s0)
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
                [mu rr wt] = approx_poles(V(:,1:k),H(1:k,1:k),h_next);
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
     function [mu rr wt tf_terms] = approx_poles(V,H,h_next)
        [Z D] = eig(H);
        ritz_vals = diag(D);
        
%         f = c'*V*Z;
%         g = Z \ (V'*r0);
%         fg = f' .* g;
%         tf_terms = bsxfun(@rdivide, fg, 1-ritz_vals*(s-s0));
%         wt = max(abs(tf_terms),[],2); %  weights
        wt = [];
        rr = abs(h_next*Z(end,:))'./abs(ritz_vals); % relative residuals
        mu = s0 + 1./ritz_vals; % poles
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
