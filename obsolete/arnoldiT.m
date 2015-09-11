function [V v_next H h_next] = arnoldiT(A,Y,r0,kmax)
% [V H h_next] = arnoldi(A,r0,kmax)
%
%  This is the Arnoldi process for computing a basis for the k-th Krylov 
%  subspace K_k(A,r0),  where k = min(d(A,r0), kmax). 
% Output:  k x k matrix H and basis vectors as columns of n x k matrix V.

tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r0);
if ~exist('kmax','var')
    kmax = n;
end

beta = norm(r0);
if abs(beta) < tol, return;  end  % no point in doing this if beta = 0;

V = zeros(n,kmax);  
H = zeros(kmax,kmax);

V(:,1) = r0/beta;
l = size(Y,2);
k = 1;
while true
   q = A(V(:,k)); 
   hk = zeros(k,1);
   
   % Orthogonalize against y's
   for i = 1:l
       yi = Y(:,i);
       q = q - (yi' * q)*yi;
   end
      
   % Orthogonalize against previous v's
   for i = 1:k
       vi = V(:,i);
       hk(i) = vi' * q;
       q = q - hk(i)*vi;
   end
   H(1:k,k) = hk;
   
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

