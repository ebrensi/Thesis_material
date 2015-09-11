function [H_hat V_hat] = RA2(A,r0,k0,num_restarts)
% [H V] = arnoldi(A,r0,k0,num_restarts)
%
%  This is the restarted Arnoldi process for computing a partial basis for the
%  k-th Krylov subspace K_k(A,r0),  where k = min(d(A,r0), kmax). 
% Output:  k x k matrix H and basis vectors as columns of n x k matrix V.

tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r0);
eta = zeros(num_restarts+1,1);

beta = norm(r0);
if abs(beta) < tol return;  end  % no point in doing this if beta = 0;

V = zeros(n, k0,num_restarts);  % this may be too much.
H = zeros(k0,k0,num_restarts);
r = 1;
v1 = r0/beta;

while r <= num_restarts
    V(:,1,r) = v1;
    k = 1;
    while true
        q = A(V(:,k,r));
        hk = zeros(k,1);
        for i = 1:k
            vi = V(:,i);
            hk(i) = q' * vi;
            q = q - hk(i)*vi;
        end
        H(1:k,k) = hk;

        h_next = norm(q);
        if h_next <= tol || k == k0
            break;      % exit the loop because k = d(A,r0)
        end
        % note that if k = kmax, we won't get here
        H(k+1,k,r) = h_next;
        V(:,k+1,r) = q / h_next;
        k = k + 1;
    end

    if k < k0
        V(:,k+1:end,r) = [];
        H(:,k+1:end,r) = [];
        H(k+1:end,:,r) = [];
    elseif size(V,2) > k0
        V(:,k0+1:end,r) = [];
    end
    
    if h_next <= tol  
        break;
    end
    v1 = q / h_next;
    r = r + 1;
    eta(r) = h_next;
end

if r < num_restarts 
    H(:,:,r+1:end) = [];
    V(:,:,r+1:end) = [];
end

V_hat = reshape(V,n,size(V,2)*size(V,3));
J_hat =  V_hat.' *  V_hat;
J_hat(abs(J_hat)<1e-10) = 0;

H_hat = zeros(r*k0);
for k = 0:r-1
    pos = k*k0+1;
    H_hat(pos:pos+k0-1,pos:pos+k0-1) = H(:,:,k+1); 
    if k < r-1
        H_hat(pos+k0,pos+k0-1) = eta(k+2);
    end
end

B_hat = J_hat*H_hat;
bo = B_hat;
qq = V_hat' * q;
ran = 1:k0*(r-1);
B_hat(ran,end) =  B_hat(ran,end) + qq(ran);
