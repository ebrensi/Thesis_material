function [V v_next] = arnoldiR(A,r0,m)
% V = arnoldiR(A,v_start,m)
%
%  This is the Arnoldi process for computing a basis for the k-th Krylov
%  subspace K_k(A,r0),  where k = min(d(A,r0), kmax).
% Output:  k x k matrix H and basis vectors as columns of n x k matrix V.

% This version of Arnoldi produces a real basis V even if A is not strictly
% real, keeping the complex vectors and using the standar
% euclidean inner product for orthogonalizing.

arnoldi_tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r0);
v_start = r0;

V = zeros(n,m*2);
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
    
    % Orthogonalize [Re(q) Im(q)] against previous v's
    q = orthagainst([real(q) imag(q)],V(:,1:kv));
    
    q1 = q(:,1);
    if norm(q1) > arnoldi_tol
        kv = kv+1;
        V(:,kv) = q1;
    end

    q2 = q(:,2);
    if norm(q2) > arnoldi_tol
        kv = kv+1;
        V(:,kv) = q2;
    end
end

if kv < m*2
    V(:,kv+1:end) = [];
end

