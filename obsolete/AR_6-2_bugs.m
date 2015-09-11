function [M H V] = AG(data,k0,num_restarts,samps)
C = []; G = []; b = []; c = []; Tfunc_real = []; L = [];

ritz_tol = 0.01;
%warning off all
if strcmp(data,'1a')
    load('example1a');
    load data1a
elseif strcmp(data,'1b')
    load('example1b');
    load data1b
else
    load(data);
end

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
eig_buffer = [];
k_tot = 0;
while r <= num_restarts
    nb = size(eig_buffer,2);
    kall = k0+nb;
    V = zeros(n, kall);  
    H = zeros(kall);
    k = 1;
    
    while true
        if k <= nb
            V(:,k) = eig_buffer(k);
        elseif k == nb+1
            V(:,k) = v_start;
        end
        q = A(V(:,k));
        hk = zeros(k,1);
        for i = 1:k
            vi = V(:,i);
            hk(i) = q' * vi;
            q = q - hk(i)*vi;
        end
        H(1:k,k) = hk;

        if k > nb
            h_next = norm(q);

            if any(k_tot == samps)
                analysis('eigs');
            end
            if h_next <= tol || k == kall
                % exit the loop because k = d(A,r0) or we've reached
                % the end of this block
                break;
            end
            % note that if k = kmax, we won't get here
            H(k+1,k) = h_next;
            V(:,k+1) = q / h_next;
            k_tot = k_tot + 1;
        end
        
        k = k + 1;

    end

    if h_next <= tol
        break;
    end
    
    v_next = q / h_next;
    eta(r) = h_next;    
    eig_buffer = check_eigs;
    v_start = v_next;
    r = r + 1;
end

r = r-1;
analysis('tfunc');

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function analysis(vis)
        Vk = V(:,1:k);
        Hk = H(1:k,1:k);
        Ik = eye(k);

        [Z D] = eig(Hk);

        br = Vk'*r0;
        cr = Vk'*c;

        g = Z.'*cr;
        ff = Z \ br;
        gf = g .* ff;
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

            weight = log(abs(gf ./ (1 - s_s0(end)*LG)));
%             weight = log(abs(gf));
%             weight = abs(LG);
            if strcmp(vis,'eigs')
                % -- Plot eigenvalues --
                scatter(real(L),imag(L),'k.')
                hold on
                ra = axis;
                scatter(real(LG),imag(LG),10,weight,'filled');

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

    function good_vecs = check_eigs
        Vk = V(:,1:k);
        last_ritz = eig(H(1:k-1,1:k-1));
        [Z D] = eig(H(1:k,1:k));
        ritz = diag(D); 
        
        jr  = zeros(2*k-1,1);
        ind = zeros(2*k-1,1);
        jr(1:k-1)  = last_ritz;
        jr(k:end)  = ritz;
        ind(k:end) = 1;
        
        [dummy idx] = sort(imag(jr));
        s_jr = jr(idx);
        d = abs(diff(s_jr)) ./ abs(s_jr(2:end));
        d_ind = diff(ind(idx));
        
%         figure;
%         plot(last_ritz,'b.')
%         hold on
%         plot(L,'k.')
%         plot(ritz,'r.')
%         hold off
        
        good = d < ritz_tol & d_ind ~= 0;
        bk = [(d_ind == -1 & good); 0];
        fd = [0; (d_ind == 1 & good)];
        goodinds = idx(bk | fd);
        goodinds = goodinds - k+1;
        goodritzvals = ritz(goodinds);
        good_ritz_vecs = Z(:,goodinds);
        
        if ~isempty(goodritzvals)
            % -- Plot eigenvalues --
            scatter(real(L),imag(L),'k.')
            hold on
            ra = axis;
            scatter(real(goodritzvals),imag(goodritzvals),10,'r','filled');
            axis(ra);
            hold off
        end
        
%         br = Vk'*r0;
%         cr = Vk'*c;
% 
%         g = Z.'*cr;
%         ff = Z \ br;
%         gf = g .* ff;   

good_vecs = Vk * good_ritz_vecs;
    end

end
