function [H V] = AT(data,m,ell,K,samps)
% [H V] = AT(data,m,ell,K,samps)

% This the first version of the full thick restart as specified 
%  by M. Eiermann et al.
%warning off all
C = []; G = []; b = []; c = []; Tfunc_real = []; L = [];

if strcmp(data,'1a')
    load('example1a');
    load data1a
elseif strcmp(data,'1b')
    load('example1b');
    load data1b
else
    load(data);
end
if isempty(samps)
    samps(1) = inf;
end
A_eigs = L;
spect_radius_A = max(abs(A_eigs));
ritz_tol = 1e-5;
%s0 = (1+sqrt(-1))*pi * 1e10 ;
s0 = pi * 1e10 ;
f = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;

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
    [V v_next H h_next] = arnoldi(Y,U,h_next,u,v_start,m);
    V_hat(:,(r-1)*m+1:r*m) = V(:,size(Y,2)+1:end);
    if h_next <= tol
        break;
    end
    
    if r < K
        [Y U X u] = deflate(V,H,h_next,ell);
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
H = CHCT*inv(CCT);
time = toc
plot_tfunc;
time
% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function [Y U X u] = deflate(V,H,h_next,ell)
        ellm = size(V,2);
        [W D] = eig(H);
        [S R] = schur(H);
        % In general, eigenvalues returned by eig and schur are not in 
        %  the same order, so we need to sort them.
        [ritz_vals Spos] = sort(ordeig(R));
        [ritz_vals Dpos] = sort(diag(D));
        err = abs(h_next*W(ellm,Dpos))'./abs(ritz_vals); % relative residuals
        bad = real(ritz_vals)>=0 | abs(ritz_vals) > spect_radius_A;
        err(bad) = inf; % we don't want these to be in the top ell values
        
        % now we sort by increasing error
        [sorted_err err_idx] = sort(err);
        figure;semilogy(sorted_err);title(sprintf('r = %d',r));
        % get positions in the schur decomp corresponding to lowest err 
        good_Spos = Spos(err_idx(1:ell));
        good = false(ellm,1);
        good(good_Spos) = true; 
        
        [S R] = ordschur(S,R,good);
        X = S(:,1:ell);

        U = R(1:ell,1:ell);
        Y = V * X;
        u = X(ellm,:)';
%         figure;
%         plot_eigs(ordeig(U),'+');
    end

    function plot_tfunc
        km = length(H);
        bk = V'*r0;
        ck = V'*c;
        J = makeJ(V,K,m,K*m);
        rank_J = rank(J)
        if rank_J < km
            warning('MATLAB:oy','J := V''*V is rank deficient (rank %d)',rank_J);
        end

        ckCCT = ck'*CCT;

        Fvec = V' * h_next*v_next;
        Fvec(m*(K-1)+1:K*m) = 0;
        %opts.POSDEF = true;
        opts.SYM = true;
        X = linsolve(J,[Fvec bk],opts);
        J_F  = X(:,1);
        %         J_F(m*(K-1)+1:K*m) = 0;
        J_bk = X(:,2);
        B = CHCT;
        B(:,km) = B(:,km) + J_F;

%         [Z D] = eig(B,CCT);
%         g = ckCCT * Z;
%         ff = Z \ J_bk;
%         gf = g' .* ff;
%         ritz_vals = diag(D);
        T1 = zeros(length(s),1);
        Tfunc = zeros(length(s),1);
        for j = 1:length(s)
            Tfunc(j)  = abs(ckCCT * ((CCT - s_s0(j)*B)\J_bk));
            %T1(j) =  abs(sum(gf ./ (1 - s_s0(j)*ritz_vals)));
        end
        %         plot it on a log scale
        %         h = figure%('visible','off');
        loglog(f,Tfunc_real,'r',f,Tfunc,f,T1)
        legend('real',...
            sprintf('orig (E=%.4g)',norm(Tfunc-Tfunc_real,inf)),...
            sprintf('eigsum (E=%.4g)',norm(T1-Tfunc_real,inf)) );
        xlabel('f','fontsize',12);
        ylabel('|H_K(s)|','fontsize',12,'Rotation',90);
        title(sprintf('m = %d, K = %d, l = %d, time: %3.4g sec',m,K,ell,time));
        drawnow
        %         saveas(h,sprintf('A1_%s_%d.png',data,k));
        %         close h
    end

    function plot_eigs(ritz_vals,marker_type)
        % -- Plot eigenvalues and ritzvals--
        scatter(real(A_eigs),imag(A_eigs),'k.')
        hold on
        ra = axis;
        scatter(real(ritz_vals),imag(ritz_vals),10,'r','filled',marker_type);
        axis(ra);
        hold off
    end

    function [V v_next H h_next] = arnoldi(Y,U,h,u,v_start,m)
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
                title(sprintf('r = %d, k = %d,',r,k-ell));
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
