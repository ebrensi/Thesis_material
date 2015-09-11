function [H V] = RA(A,r0,k0,restarts)
% [H V] = arnoldi(A,r0,k0,restarts)
%
%  This is the restarted Arnoldi process for computing a partial basis for the
%  k-th Krylov subspace K_k(A,r0),  where k = min(d(A,r0), kmax). 
% Output:  k x k matrix H and basis vectors as columns of n x k matrix V.

tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r0);

beta = norm(r0);
if abs(beta) < tol return;  end  % no point in doing this if beta = 0;

V = zeros(n, k0,restarts);  % this may be too much.
H = zeros(k0,k0,restarts);
r = 1;
v1 = r0/beta;

while r <= restarts
    V(:,1,r) = v1;
    k = 1;
    while true
        Vk = V(:,1:k,r);
        vk = Vk(:,k);
        q = A(vk);

        hk = Vk' * q;
        H(1:k,k,r) = hk;
        q = q - Vk * hk;

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
end

if r < restarts 
    H(:,:,r+1:end) = [];
    V(:,:,r+1:end) = [];
end
