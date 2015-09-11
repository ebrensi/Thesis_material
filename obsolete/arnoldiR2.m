function [V V2 Vc] = arnoldiR2(A,r,m)
% V = arnoldiR(A,v_start,m)
%
%  This is the real-Arnoldi process for computing a basis for the k-th Split-Krylov
%  subspace K_k(A,r)*.

% This version of arnoldiR uses standard complex Arnoldi to compute a basis
%  for K_k(A,r), splits the basis into real and imaginary parts, and 
%  orthonormalizes them.

arnoldi_tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r);
v_start = r;

% V = zeros(n,2*m);
V2 = zeros(n,2*m);
Vc = zeros(n,m);

k = 0;
kv = 0;
while  k < m
    if k
        q = A(v_next);
    else
        q = v_start;
    end
    
    % Orthogonalize q against previous vc's
    v_next = orthagainst(q,Vc(:,1:k));
    k = k + 1;
    Vc(:,k) = v_next;
    
    V2(:,2*k-1:2*k) = [real(v_next) imag(v_next)];
%     V2(:,2*k-1:2*k) = [real(q) imag(q)];
    
    % Orthogonalize [Re(q) Im(q)] against previous v's
%     q = orthagainst([real(q) imag(q)],V(:,1:kv));
%     
%     q1 = q(:,1);
%     if norm(q1) > arnoldi_tol
%         kv = kv+1;
%         V(:,kv) = q1;
%     end
% 
%     q2 = q(:,2);
%     if norm(q2) > arnoldi_tol
%         kv = kv+1;
%         V(:,kv) = q2;
%     end
    
end

% if kv < m*2
%     V(:,kv+1:end) = [];
% end

V = orthagainst(V2);

