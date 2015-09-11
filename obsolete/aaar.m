function aaar(data,k0,num_restarts)
C = []; G = []; b = []; c = []; Tfunc_real = [];
load(data);
s0 = pi * 1e10;
f = logspace(8,10,200);  % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;
% Tfunc_real = zeros(length(s),1);
% for j = 1:length(s)
% 	Tfunc_real(j) = abs(c.' * ((s(j)*G - C)\b));
% end
load data1a 
plot(L,'r.')
C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);
eta = zeros(num_restarts+1,1);

beta = norm(r0);
if abs(beta) < tol return;  end  % no point in doing this if beta = 0;

V = zeros(n, k0,num_restarts);  % this may be too much.
H = zeros(k0,k0,num_restarts);
r = 1;
v1 = r0/beta;
kk = 0;
while r <= num_restarts
    V(:,1,r) = v1;
    kk = kk + 1;
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
        if h_next <= tol || k == k0
            break;      % exit the loop because k = d(A,r0)
        end
        % note that if k = kmax, we won't get here
        H(k+1,k,r) = h_next;
        V(:,k+1,r) = q / h_next;
        k = k + 1;
        kk = kk + 1;
    end
    
    if h_next <= tol  
        break;
    end
    
    if any(kk == [50 100])
        Hrplot;
    end
    
    v1 = q / h_next;
    r = r + 1;
    eta(r) = h_next;
end


% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function Hrplot
        rk0 = r*k0;
        Vr = reshape(V(:,:,1:r),n,rk0);
        Jr =  Vr.' *  Vr;
        Jr(abs(Jr)<eps) = 0;

        Hr = zeros(rk0);
        for k = 0:r-1
            pos = k*k0+1;
            Hr(pos:pos+k0-1,pos:pos+k0-1) = H(:,:,k+1);
            if k < r-1
                Hr(pos+k0,pos+k0-1) = eta(k+2);
            end
        end
        %Hr(abs(Hr)<eps) = 0;
        Br = Jr*Hr;
        qq = Vr' * q;
        ran = 1:k0*(r-1);
        Br(ran,end) =  Br(ran,end) + qq(ran);

        e1 = eye(rk0,1);
        br = Vr'*r0;
        cr = Vr'*c;
        
		%% - alternatively -
		Gr = Hr;
		Gr(:,end) = Gr(:,end) + Jr\qq;
		LG = eig(Gr);
		figure;
		plot(LG,'.');
		%%
		
        Tfunc = zeros(length(s),1);
		Tfunc = zeros(length(s),1);
		Ir = eye(r*k0); 
		for j = 1:length(s)
            Tfunc(j)  = abs(cr' * ((Jr - (s(j)-s0)*Br)\br));
			Tfunc2(j) = abs(cr' * ((Jr*(Ir - (s(j)-s0)*Gr))\br));
        end
%         plot it on a log scale
%         h = figure('visible','off');
        h = figure;
        if exist('Tfunc_real','var')
            loglog(f,Tfunc_real,'r',f,Tfunc,f,Tfunc2)
			legend('real','orig','alt')
        else
            loglog(f,Tfunc);
        end

        xlabel('f','fontsize',12);
        ylabel('|Hr(s)|','fontsize',12,'Rotation',90');
        title(sprintf('r = %d, k = %d',r,r*k0));
%         saveas(h,sprintf('A1_%s_%d.png',data,k));
%         close h


    end
end
