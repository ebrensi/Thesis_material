function [H V] = AT(data,m,min_ell_ess,K,samps)
% [H V] = AT(data,m,ell,K,samps)

% This uses thick restart with variable l.  
% version 2 used prior to 9/30/2008
% this version abandoned because it doesn't deal with complex s0 properly.
C = []; G = []; b = []; c = []; Tfunc_real = []; L = [];

if strcmp(data,'1a')
    load('ex1a');
    load data1a
    A_eigs = L;
elseif strcmp(data,'1b')
    load('ex1b');
    load data1b
    A_eigs = L;
else
    clear Tfunc_real L
    load(data);
end
if isempty(samps)
    samps(1) = inf;
end
ritz_tol = 1e-5;
s0 = (1+sqrt(-1))*pi * 1e10 ;
% s0 = pi * 1e10 ;
frq= logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;

s_s0 = s - s0;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);

sc = 1;
k_tot = 0;
v_start = r0/norm(r0);
Y = []; U = []; u = 1; h_next = 1;
CCT = zeros(K*m);
CHCT = zeros(K*m);
C = eye(m);
V_hat = zeros(n,K*m);

r = 1;
tic;
while r <= K 
    [V v_next H h_next] = thickstart_arnoldi(Y,U,h_next,u,v_start,m);
    V_hat(:,(r-1)*m+1:r*m) = V(:,size(Y,2)+1:end);
    if h_next <= tol
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

if h_next <= tol
   fprintf('The Krylov space is maxed out.\n\n'); 
end
r = r-1;

V = V_hat;
time = toc;
plot_tfunc;
% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function [Y U X u ell] = deflate(V,H,h_next,min_ell_ess)
        ellm = size(V,2);
        [Z D] = eig(H);
        [S R] = schur(H);
        % In general, eigenvalues returned by eig and schur are not in 
        %  the same order, so we need to sort them.
        [ritz_vals Spos] = sort(ordeig(R));
        [ritz_vals Dpos] = sort(diag(D));
        Z = Z(:,Dpos);
        ess = imag(ritz_vals) <= 0;
        ess_ritz = ritz_vals(ess);
        cplx = imag(ess_ritz) < 0;
        
        % ------ 1st criterion (relative residual error) -----
        rr = abs(h_next*Z(end,ess))'./abs(ess_ritz); % relative residuals
        [sorted_rr idx] = sort(rr); st = 'residual';
        % -----------------------------------------------------

        
        % ------ 2nd criterion (relevance of term in tfunc sum) -----
        g = Z'*V'*c;
        f = Z \ (V'*r0);
        gf = conj(g) .* f; % g must be conjugated. recall H(s) = (g'(I-(s-s0)A)^-1)f)
        tf_terms_all = bsxfun(@rdivide, gf, (1 - ritz_vals*s_s0));
        tf_terms = tf_terms_all(ess,:);
        tf_terms(cplx,:) = tf_terms(cplx,:) + tf_terms_all(~ess,:);
        clear tf_terms_all;
        wt = max(abs(tf_terms),[],2);
        [sorted_wt idx] = sort(wt,'descend'); st = 'weight'
        % -----------------------------------------------------------
        
        [score idx] = sort(wt./rr,'descend'); st = 'score';

         good_ess = rr < ritz_tol;
         if nnz(good_ess) < min_ell_ess
             good_ess(idx(1:min_ell_ess)) = true;
         end

        good = false(ellm,1);        
        good(ess) = good_ess;
        good(~ess) = good_ess(cplx);
        ell = nnz(good)
        
        [S R] = ordschur(S,R,good);
        X = S(:,1:ell);
        U = R(1:ell,1:ell);
        Y = V * X;
        u = X(end,:)';
       
%         figure;
%         plot_eigs(ritz_vals(good),'+');
    end

    function plot_tfunc
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
              
        Tfunc = zeros(length(s),1);

        for j = 1:length(s)
            Tfunc(j)  = abs(ckCCT * ((CCT - s_s0(j)*B) \ J_bk));
        end
        %         plot it on a log scale
        %         h = figure%('visible','off');
        if exist('Tfunc_real','var')
            loglog(frq,Tfunc_real,'r',frq,Tfunc);
            err = norm(Tfunc-Tfunc_real,inf) / norm(Tfunc_real,inf);
            legend('real', sprintf('orig (rel err=%.4g)',err) );
        else
            loglog(frq,Tfunc);
        end
        xlabel('frq','fontsize',12);
        ylabel('|H_K(s)|','fontsize',12,'Rotation',90);
        title(sprintf('m = %d, K = %d, l_{min} = %d, time: %3.4g sec',m,K,min_ell_ess,time));
        drawnow
        %         saveas(h,sprintf('A1_%s_%d.png',data,k));
        %         close h
    end

    function plot_eigs(ritz_vals,marker_type, weights)
        % -- Plot eigenvalues and ritzvals--
        if exist('A_eigs','var');
            scatter(real(A_eigs),imag(A_eigs),'k.');
        end
        hold on
        ra = axis;
        if nargin == 1
            scatter(real(ritz_vals),imag(ritz_vals),10,'r','filled');
        elseif nargin == 2
            scatter(real(ritz_vals),imag(ritz_vals),10,'r',marker_type);
        elseif nargin == 3
            scatter(real(ritz_vals),imag(ritz_vals),10,weights,'filled',marker_type);
        end
        axis(ra);
        hold off
    end

    function [V v_next H h_next] = thickstart_arnoldi(Y,U,h,u,v_start,m)
        l = size(Y,2);
        k =l+1;
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
        while  k <= ellm && h_next > tol
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
            
            % **** visualization of approxmate e-values
            k_tot = k_tot + 1;
            if k_tot == samps(sc)
                plot_eigs(eig(H(1:k,1:k)),'o');
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
