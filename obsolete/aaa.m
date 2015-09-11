function aaa(data,kmax)
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

s0 = pi * 1e10;
f = logspace(8,10,200);  % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*f;
% Tfunc_real = zeros(length(s),1);
% for j = 1:length(s)
% 	Tfunc_real(j) = abs(c.' * ((s(j)*G - C)\b));
% end

C_s0G = C-s0*G;
r0 = -C_s0G\b;

tol = eps;
n = length(r0);
beta = norm(r0);
if abs(beta) < tol, return;  end  % no point in doing this if beta = 0;

V = zeros(n,kmax);  
H = zeros(kmax,kmax);

V(:,1) = r0/beta;
k = 1;
while true
    q = A(V(:,k));
    hk = zeros(k,1);
    for i = 1:k
       vi = V(:,i);
       hk(i) = q' * vi;
       q = q - hk(i)*vi; 
    end
    H(1:k,k) = hk;

    h_next = norm(q);
    if h_next <= tol || k > kmax
        break;      % exit the loop because k = d(A,r0)
    end
    % note that if k = kmax, we won't get here
    H(k+1,k) = h_next;
    V(:,k+1) = q / h_next;
    k = k + 1;

    if any(k == [50 60 70 80 90 100])
		LH = eig(H(1:k,1:k));
		figure;
% 		plot(LH,'.');
        Hkplot(H(1:k,1:k),V(:,1:k),k)
    end
end

if k < kmax
    V(:,k+1:end) = [];
    H(:,k+1:end) = [];
    H(k+1:end,:) = [];
elseif size(V,2) > kmax
    V(:,kmax+1:end) = [];
end
figure;
plot(L,'r.');

% -------------------------------------------------------------
    function q = A(v)
        q = C_s0G\(G*v);
    end

    function Hkplot(Hk,Vk,k)
        e1 = [1; zeros(k-1,1)];
        bk = norm(r0)*e1;
        ck = Vk'*c;
        Ik = eye(k);
        Tfunc = zeros(length(s),1);
        for j = 1:length(s)
            Tfunc(j) = abs(ck' * ((Ik - (s(j)-s0)*Hk)\bk));
        end
%         plot it on a log scale
%         h = figure('visible','off');
        h = figure;
        if exist('Tfunc_real','var')
            loglog(f,Tfunc_real,'r',f,Tfunc)
        else
            loglog(f,Tfunc);
        end

        xlabel('f','fontsize',12);
        ylabel('|H_k(s)|','fontsize',12,'Rotation',90');
        title(sprintf('k = %d',k));
%         saveas(h,sprintf('A1_%s_%d.png',data,k));
%         close h
    end
end
