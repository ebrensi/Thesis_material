function [M H V] = AG(data,k0,num_restarts,samps)
C = []; G = []; b = []; c = []; Tfunc_real = []; L = [];
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

V = zeros(n, k0,num_restarts);  % this may be too much.
H = zeros(k0,k0,num_restarts);

beta = norm(r0);
if abs(beta) < tol
    return;
end  % no point in doing this if beta = 0;

% set(gca,'nextplot','replacechildren'); % set up animation
frm = 1;
r = 1;
v1 = r0/beta;
k_tot = 0;
while r <= num_restarts
    V(:,1,r) = v1;
    k_tot = k_tot + 1;
    k = 1;
    while true
        q = A(V(:,k,r));
        hk = zeros(k,1);
        for i = 1:k
            vi = V(:,i,r);
            hk(i) = q' * vi;
            q = q - hk(i)*vi;
        end
        H(1:k,k,r) = hk;

        h_next = norm(q);

        if any(k_tot == samps)
            analysis('poles');
        end
        if h_next <= tol || k == k0
            % exit the loop because k = d(A,r0) or we've reached
            % the end of this block
            break;
        end
        % note that if k = kmax, we won't get here
        H(k+1,k,r) = h_next;
        V(:,k+1,r) = q / h_next;

        k = k + 1;
        k_tot = k_tot + 1;
    end

    if h_next <= tol
        break;
    end

    v1 = q / h_next;
    r = r + 1;
    eta(r) = h_next;
end

r = r-1;
analysis('tfunc');

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function analysis(vis)
        persistent lastLG
        Vr = [reshape(V(:,:,1:r-1),n,(r-1)*k0) V(:,1:k,r)];
        Jr =  makeJr(Vr,r,k0,k_tot);
        %[k_tot cond(Jr)]
        Hr = makeHr(H,r,k0,k_tot);

        Br = Jr*Hr;
        qq = Vr.' * q;
        zero_ran = k0*(r-1)+1:k_tot;
        qq(zero_ran) = 0;
        Br(:,end) =  Br(:,end) + qq;
        [Z D] = eig(Br,Jr);

        br = Vr'*r0;
        cr = Vr'*c;
        % -----------------
        g = Z.'*cr;
        ff = (Jr*Z) \ br;
        gf = g .* ff;
        LG = diag(D); 
                    
%             S = s0 + 1./LG;
%             ritz_mask = smin < abs(S) & abs(S) < smax;
%             LG = LG(ritz_mask);
%             gf = gf(ritz_mask);
            % -----------------
        if strcmp(vis,'tfunc')
            I = eye(k_tot);
            Tfunc = zeros(length(s),1);
            T1 = zeros(length(s),1);
            for j = 1:length(s)
                Tfunc(j)  = abs(cr' * ((Jr - s_s0(j)*Br)\br));
                T1(j) =  abs(sum(gf ./ (1 - s_s0(j)*LG)));
                % Tfunc2(j) = abs(cr' * ((Jr*(Ir - s_s0(j)*Gr))\br));
            end
            %         plot it on a log scale
            %         h = figure('visible','off');

            h = figure;
            if exist('Tfunc_real','var')
                loglog(f,Tfunc_real,'r',f,Tfunc,f,T1)
                %legend('real', sprintf('orig (E=%.4g)',norm(Tfunc-Tfunc_real,inf)));
                legend('real', sprintf('orig (E=%.4g)',norm(Tfunc-Tfunc_real,inf)), sprintf('alt (E=%.4g)',norm(T1-Tfunc_real,inf)) );
            else
                loglog(f,Tfunc);
            end

            xlabel('f','fontsize',12);
            ylabel('|Hr(s)|','fontsize',12,'Rotation',90');
            title(sprintf('r = %d, k_{tot} = %d,',r,k_tot));
            %         saveas(h,sprintf('A1_%s_%d.png',data,k));
            %         close h
        end



        if strcmp(vis,'eigs') || strcmp(vis,'poles')
            LG = diag(D);
            
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
        lastLG = LG;
    end

    function Jr = makeJr(V,r,k0,k_tot)
        Jr = eye(k_tot);
        if r > 1
            for j = 1:r-1
                block_j = (j-1)*k0+1:j*k0;
                after_j = j*k0+1:k_tot;
                Jr(block_j,after_j) = V(:,block_j)' * V(:,after_j);
            end
            Jr = Jr + triu(Jr,1).';
        end
    end


    function  Hr = makeHr(H,r,k0,k_tot)
        Hr = zeros(k_tot);
        if r > 1
            for j = 1:r-1
                block_j = (j-1)*k0+1:j*k0;
                Hr(block_j,block_j) = H(:,:,j);
                Hr(block_j(end)+1,block_j(end)) = eta(j+1);
            end
        end
        block_r = (r-1)*k0+1:k_tot;
        ki = length(block_r);
        Hr(block_r,block_r) = H(1:ki,1:ki,r);
    end

end
