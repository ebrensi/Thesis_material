function [H V] = ARS(data,m,ell,num_runs,samps)
% [M H V] = ARS(data,m,ell,num_runs,samps)
% This is the transfer function approximation using thick restarting 
%  but not keeping V and H from previous restarts. 7/17/2008
%  another difference between this and previous versions is that 
%  we keep the best ell ritz values regardless of convergence.
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
%s0 = (1+sqrt(-1))*pi * 1e10 ;
s0 = pi * 1e10 ;
f = logspace(8,10,200);  % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;

s_s0 = s - s0;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);

sc = 1;
h_next = norm(r0);
v_start = r0/h_next;

rstrt_k = m;
V = zeros(n,m+1);
H = zeros(m+1);
k_tot = 1;
k = 1;
r = 1;
while r <= num_runs
    V(:,k) = v_start;
    while  h_next > tol && k <= rstrt_k
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
        
        H(k+1,k) = h_next;
        V(:,k+1) = q / h_next;
       
        if k_tot == samps(sc)
            plot_eigs(eig(H(1:k,1:k)),'o');
            title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
            drawnow;
            if sc < length(samps)
                sc = sc+1;
            end
        end
        
        k_tot = k_tot + 1;
        k = k + 1;      
    end

    if h_next <= tol
        break;
    end
    
    k = rstrt_k;
    if r < num_runs
        v_start = V(:,rstrt_k+1);
        
        [V_hat H_hat s_vec] = deflate(V,H,rstrt_k,ell);
        ell = length(H_hat);
        rstrt_k = m+ell;
        V = zeros(n, rstrt_k+1);
        H = zeros(rstrt_k+1);
        if ell
            V(:,1:ell) = V_hat;
            H(1:ell,1:ell) = H_hat; 
            H(ell+1,1:ell) = h_next*s_vec';
        end
        k = ell + 1;
    end
    r = r + 1;
end

r = r-1;
k_tot = k_tot-1;
plot_tfunc(k);

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function [V_hat H_hat s_vec] = deflate(V,H,m,ell)
        Hm = H(1:m,1:m); 
        Vm = V(:,1:m);
        [W D] = eig(Hm);
        [U R] = schur(Hm);
        % In general, eigenvalues returned by eig and schur are not in 
        %  the same order, so we need to sort them.
        [sorted_ritz Spos] = sort(ordeig(R));
        [sorted_ritz Dpos] = sort(diag(D));
        err = abs(H(m+1,m)*W(m,Dpos))'./abs(sorted_ritz); % relative residuals
        bad = real(sorted_ritz)>=0 | abs(sorted_ritz) > spect_radius_A;
        err(bad) = inf; % we don't want these to be in the top ell values
        
        % now we sort by increasing error
        [sorted_err err_idx] = sort(err);
        % get positions in the schur decomp corresponding to lowest err 
        good_Spos = Spos(err_idx(1:ell));
        good = false(m,1);
        good(good_Spos) = true; 
              
        [U R] = ordschur(U,R,good);
        U_hat = U(:,1:ell);
        H_hat = R(1:ell,1:ell);
        V_hat = Vm * U_hat;
        s_vec = U_hat(m,:)';

%         figure;
%         plot_eigs(ordeig(H_hat),'+');
    end

    function plot_tfunc(k)
        Vk = V(:,1:k);
        Hk = H(1:k,1:k);
        Ik = eye(k);
        br = Vk'*r0;
        cr = Vk'*c;      
        
        [Z D] = eig(Hk);
        g = Z.'*cr;
        ff = Z \ br;
        gf = g .* ff;
        ritz_vals = diag(D);
        Tfunc = zeros(length(s),1);
        T1 = zeros(length(s),1);
        for j = 1:length(s)
            Tfunc(j)  = abs(cr' * ((Ik - s_s0(j)*Hk)\br));
            T1(j) =  abs(sum(gf ./ (1 - s_s0(j)*ritz_vals)));
        end
        %         plot it on a log scale
        %         h = figure%('visible','off');
        loglog(f,Tfunc_real,'r',f,Tfunc,f,T1)
        legend('real',...
            sprintf('orig (E=%.4g)',norm(Tfunc-Tfunc_real,inf)),...
            sprintf('eigsum (E=%.4g)',norm(T1-Tfunc_real,inf)) );
        xlabel('f','fontsize',12);
        ylabel('|H_k(s)|','fontsize',12,'Rotation',90);
        title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
        %         saveas(h,sprintf('A1_%s_%d.png',data,k));
        %         close h
        
    end

    function plot_eigs(ritz_vals,marker_type)
        % -- Plot eigenvalues and ritzvals--
        scatter(real(A_eigs),imag(A_eigs),'k.')
        hold on
        ra = axis;
        scatter(real(ritz_vals),imag(ritz_vals),15,'r','filled',marker_type);
        axis(ra);
        hold off
    end

end
