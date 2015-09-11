function [H V] = ARS(data,m,num_runs,samps)
% [M H V] = AR(data,m,num_runs,samps)
% This is the first version of AR that uses shur decompostion 
%  for thick restart vectors. 
%  This version is getting ditched because I've decided to use only 
%    residue ritz-value convergence criterion.

C = []; G = []; b = []; c = []; Tfunc_real = []; L = [];

ritz_tol = 0.001;

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
spect_radius_A = max(abs(L));
%s0 = (1+sqrt(-1))*pi * 1e10 ;
s0 = pi * 1e10 ;
f = logspace(8,10,200);  % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;

s_s0 = s - s0;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);

% set(gca,'nextplot','replacechildren'); % set up animation
sc = 1;
frm = 1;
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
            analysis('eigs');
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
        
        [V_hat H_hat] = deflate(V,H,rstrt_k);
        ell = length(H_hat);
        rstrt_k = m+ell;
        V = zeros(n, rstrt_k+1);
        H = zeros(rstrt_k+1);
        if ell
            V(:,1:ell) = V_hat;
            H(1:ell,1:ell) = H_hat;
        end
        k = ell + 1;
    end
    r = r + 1;
end

r = r-1;
k_tot = k_tot-1;
analysis('tfunc');

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function [V_hat H_hat] = deflate(V,H,m)
        [U,R] = schur(H(1:m,1:m),'complex');
        H_eigs = ordeig(R);
        good = real(H_eigs)<=0 & abs(H_eigs)<=spect_radius_A;
        n_good = nnz(good);
        [U,R] = ordschur(U,R,good); 
        H_eigs = ordeig(R);
        
        
        % (1) ritz value convergence criterion based on change over iterations 
        %  ********************************************************
        prev_H_eigs = eig(H(1:m-1,1:m-1));
        discard_mask = real(prev_H_eigs)>0 | abs(prev_H_eigs)> spect_radius_A ;
        prev_H_eigs(discard_mask) = [];
        
        combined_eigs = [H_eigs(1:n_good); prev_H_eigs(:)];
        n_combined_eigs = length(combined_eigs);
        tag = zeros(n_combined_eigs,1);
        tag(1:n_good) = 1;
        [dummy idx] = sort(imag(combined_eigs));
        s_combined_eigs = combined_eigs(idx);
        tag_diff = diff(tag(idx));
        e_num = diff(s_combined_eigs);
        e_denom = zeros(n_combined_eigs-1,1);
        orig_idx = zeros(n_combined_eigs,1);
        mask = (tag_diff < 0);
        e_denom(mask) = s_combined_eigs(mask);
        orig_idx(mask) = idx(mask);
        mask = (tag_diff > 0);
        shift_mask = [false; mask];
        e_denom(mask) = s_combined_eigs(shift_mask);
        orig_idx(mask) = idx(shift_mask);
        e = abs(e_num) ./ abs(e_denom);
        good = e < ritz_tol;
        converged = orig_idx(good);
        % *************************************************************
        
        good = false(length(H_eigs),1);
        good(converged) = true;
        n_good = nnz(good);
        [U,R] = ordschur(U,R,good);
        
        U_hat = U(:,1:n_good);
        H_hat = R(1:n_good,1:n_good);
        V_hat = V(:,1:m) * U_hat;
        
        % alternatively,
%         [W D] = eig(H(1:m,1:m));
%         err = abs(H(m+1,m)*W(m,:))' ./ abs(ordeig(D));
%         good = err < 0.002;
%         n_good = nnz(good);
%         [U R] = ordschur(W,D,good);
%         U_hat2 = U(:,1:n_good);
%         H_hat2 = R(1:n_good,1:n_good);
        % 
        
%         [Z D] = eig(H_hat);
%         err = abs(H(m+1,m)*U_hat(m,:)*Z)' ./ abs(ordeig(D)); 
               
        if ~isempty(U_hat)
            good_ritz_vals = diag(H_hat);
            % -- Plot eigenvalues --
            figure
            scatter(real(L),imag(L),'k.')
            hold on
            ra = axis;
            scatter(real(good_ritz_vals),imag(good_ritz_vals),15,'r+');
            axis(ra);
            hold off
            figure
        end
    end

    function analysis(vis)
        Vk = V(:,1:k);
        Hk = H(1:k,1:k);
        Ik = eye(k);
        br = Vk'*r0;
        cr = Vk'*c;
        
        [Z D] = eig(Hk);

        g = Z.'*cr;
        ff = Z \ br;

        gf = g .* ff;
        LG = diag(D);

       
        if strcmp(vis,'tfunc')
            Tfunc = zeros(length(s),1);
            T1 = zeros(length(s),1);
            for j = 1:length(s)
                Tfunc(j)  = abs(cr' * ((Ik - s_s0(j)*Hk)\br));
                T1(j) =  abs(sum(gf ./ (1 - s_s0(j)*LG)));
            end
            %         plot it on a log scale
            %         h = figure('visible','off');

            h = figure;
            if exist('Tfunc_real','var')
                loglog(f,Tfunc_real,'r',f,Tfunc,f,T1)
                legend('real', sprintf('orig (E=%.4g)',norm(Tfunc-Tfunc_real,inf)), sprintf('eigsum (E=%.4g)',norm(T1-Tfunc_real,inf)) );
            else
                loglog(f,Tfunc);
            end

            xlabel('f','fontsize',12);
            ylabel('|H_k(s)|','fontsize',12,'Rotation',90');
            title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
            %         saveas(h,sprintf('A1_%s_%d.png',data,k));
            %         close h
        end



        if strcmp(vis,'eigs') 
            %weight = log(abs(gf ./ (1 - s_s0(end)*LG)));
            %             weight = log(abs(gf));
            %             weight = abs(LG);
            if strcmp(vis,'eigs')
                % -- Plot eigenvalues --
                scatter(real(L),imag(L),'k.')
                hold on
                ra = axis;
                %scatter(real(LG),imag(LG),10,weight,'filled');
                scatter(real(LG),imag(LG),10,'r','filled');
            end
            title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
            axis(ra);
            hold off
            %M(frm) = getframe(gcf);
            drawnow
            frm = frm + 1;
        end
    end

end
