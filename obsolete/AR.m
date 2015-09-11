function [M H V] = AR(data,k0,num_restarts,samps)
% [M H V] = AR(data,k0,num_restarts,samps)
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
spect_radius_A = max(abs(L));
%s0 = (1+sqrt(-1))*pi * 1e10 ;
s0 = pi * 1e10 ;
f = logspace(8,10,200);  % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;

% -- we only care about a certain frequency range? --
smin = min(abs(s));
smax = max(abs(s));
S_true = s0 + 1./L;
true_mask = smin < abs(S_true) & abs(S_true) < smax;
S_true = S_true(true_mask);
% L = L(true_mask);
% --

s_s0 = s - s0;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);
eta = zeros(num_restarts+1,1);


beta = norm(r0);
if abs(beta) < tol
    return;
end  % no point in doing this if beta = 0;

% set(gca,'nextplot','replacechildren'); % set up animation
frm = 1;
r = 1;
v_start = r0/beta;
good_vecs = [];
Y = [];
k_tot = 0;
while r <= num_restarts
    nY = size(Y,2);
    V = zeros(n, k0);
    H = zeros(k0);
    k = 1;
    V(:,1) = v_start;
    while true
        q = A(V(:,k));
        hk = zeros(k,1);
        % Orthogonalize against approximate eigenvectors 
        for i = 1:nY
            q = q - (vi' * q)*Y(:,i);
        end
        % Orthogonalize against previous v's 
        for i = 1:k
            vi = V(:,i);
            hk(i) = vi' * q;
            q = q - hk(i)*vi;
        end
        H(1:k,k) = hk;

        h_next = norm(q);

        if any(k_tot == samps)
            analysis('eigs');
        end
        if h_next <= tol || k == k0
            % exit the loop because k = d(A,r0) or we've reached
            % the end of this block
            break;
        end
        % note that if k = kmax, we won't get here
        H(k+1,k) = h_next;
        V(:,k+1) = q / h_next;
        k_tot = k_tot + 1;

        k = k + 1;
    end

    if h_next <= tol
        break;
    end

    v_next = q / h_next;
    eta(r) = h_next;
    good_vecs = [good_vecs check_eigs];
    %Y = orth(good_vecs);
    v_start = v_next;
    r = r + 1;
end

r = r-1;
analysis('tfunc');

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function good_vecs = check_eigs
        Vk = V(:,1:k);
        last_ritz = eig(H(1:k-1,1:k-1));
        discard_mask = imag(last_ritz)<0 | real(last_ritz)>0 | abs(last_ritz)> spect_radius_A ;
        last_ritz(discard_mask) = []; % keep only positive half of conjugate pairs
        n_last_ritz = length(last_ritz);

        [Z D] = eig(H(1:k,1:k));
        ritz = diag(D);
        discard_mask = imag(ritz)<0 | real(ritz)>0 | abs(ritz)> spect_radius_A;
        ritz(discard_mask) = []; % keep only positive half of conjugate pairs
        Z(:,discard_mask) = [];
        n_ritz = length(ritz);

        combined_ritz = [last_ritz(:); ritz(:)];
        n_combined_ritz = length(combined_ritz);
        tag = zeros(n_combined_ritz,1);
        tag(n_last_ritz+1:end) = 1;
        [dummy idx] = sort(imag(combined_ritz));
        s_combined_ritz = combined_ritz(idx);

        tag_diff = diff(tag(idx));
        e_num = diff(s_combined_ritz);
        e_denom = zeros(n_combined_ritz-1,1);
        orig_idx = zeros(n_combined_ritz,1);
        
        mask = (tag_diff < 0);
        e_denom(mask) = s_combined_ritz(mask);
        orig_idx(mask) = idx(mask);
        
        mask = (tag_diff > 0);
        shift_mask = [false; mask];
        e_denom(mask) = s_combined_ritz(shift_mask);
        orig_idx(mask) = idx(shift_mask);
        
        e = abs(e_num) ./ abs(e_denom);

        %         figure;
        %         plot(last_ritz,'b.')
        %         hold on
        %         plot(L,'k.')
        %         plot(ritz,'r.')
        %         hold off

        good = e < ritz_tol;
        goodinds = orig_idx(good) - n_last_ritz;
        good_ritz_vecs = Z(:,goodinds);
        cmplxvecs = any(imag(good_ritz_vecs));
        good_ritz_vecs = [good_ritz_vecs conj(good_ritz_vecs(:,cmplxvecs))];

        if ~isempty(good_ritz_vecs)
            goodritzvals = [ritz(goodinds); conj(ritz(goodinds))];
            % -- Plot eigenvalues --
            figure
            scatter(real(L),imag(L),'k.')
            hold on
            ra = axis;
            scatter(real(goodritzvals),imag(goodritzvals),15,'r+');
            axis(ra);
            hold off
        end

        W = orth(good_ritz_vecs);
        % W = good_ritz_vecs;
        good_vecs = Vk * W;
    end

    function analysis(vis)
        Vk = V(:,1:k);
        Hk = H(1:k,1:k);
        Ik = eye(k);

        [Z D] = eig(Hk);

%         br = beta*eye(k,1);
%         cr = Vk'*c;
% 
%         g = Z.'*cr;
%         ff = Z \ br;
% 
%         gf = g .* ff;
        LG = diag(D);

        %             S = s0 + 1./LG;
        %             ritz_mask = smin < abs(S) & abs(S) < smax;
        %             LG = LG(ritz_mask);
        %             gf = gf(ritz_mask);
        % -----------------
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



        if strcmp(vis,'eigs') || strcmp(vis,'poles')
            S = s0 + 1./LG;
            ritz_mask = smin < abs(S) & abs(S) < smax;
            %LG = LG(ritz_mask);
            S = S(ritz_mask);

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
            elseif strcmp(vis,'poles')
                %  -- Plot Poles --
                scatter(real(S_true),imag(S_true),'k.')
                hold on
                ra = axis;
                scatter(real(S),imag(S),10,weight(ritz_mask),'filled');
            end
            title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
            axis(ra);
            hold off
            M(frm) = getframe(gcf);
            frm = frm + 1;
        end
    end

end
