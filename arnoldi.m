function [V v_next H op_norm_est] = arnoldi(A,r0,kmax)
% [V v_next H] = arnoldi(A,r0,kmax)
%
%  This is the Arnoldi process for computing a basis for the k-th Krylov 
%  subspace K_k(A,r0),  where k = min(d(A,r0), kmax). 
% Output:  k x k matrix H and basis vectors as columns of n x k matrix V.

reorthogonalize = false;
op_norm_est = 0;
tol = eps;

if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r0);
if ~exist('kmax','var')
    kmax = n;
end

V = zeros(n,kmax);  
H = zeros(kmax,kmax);

V(:,1) = r0/norm(r0);
k = 1;
while true
   q = A(V(:,k));
   op_norm_est = max(norm(q),op_norm_est);
   hk = zeros(k,1);
   % Orthogonalize against previous v's
   for i = 1:k
       vi = V(:,i);
       hk(i) = vi' * q;
       q = q - hk(i)*vi;
   end
   H(1:k,k) = hk;
   % optional re-orthogonalization
   if reorthogonalize
	   for i = 1:k
		   vi = V(:,i);
		   h2k(i) = vi' * q;
		   q = q - h2k(i)*vi;
	   end
	   H(1:k,k) = H(1:k,k) + h2k;
   end
   
   h_next = norm(q);
   v_next = q / h_next;
   if h_next <= tol || k == kmax
       break;      % exit the loop because k = d(A,r0)
   end 
   % note that if k = kmax, we won't get here
   H(k+1,k) = h_next;
   V(:,k+1) = v_next;
   k = k + 1;
end

if k < kmax
    V(:,k+1:end) = [];
    H(:,k+1:end) = [];
    H(k+1:end,:) = [];
elseif size(V,2) > kmax
    V(:,kmax+1:end) = [];
end

