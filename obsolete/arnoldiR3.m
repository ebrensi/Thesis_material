function [V V2 Vr] = arnoldiR3(A,r,m)
% V = arnoldiR(A,v_start,m)
%
%  This is the real-Arnoldi process for computing a basis for the k-th Split-Krylov
%  subspace K_k(A,r)*.
%
% This version of arnoldiR uses equivalent real Arnoldi to compute a basis
%  for K_k(A^,r^), splits the basis, and orthonormalizes the halves.

arnoldi_tol = eps;
if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r);
v_start = r;

% V = zeros(n,2*m);
V2 = zeros(n,2*m);
Vr = zeros(2*n,m);

k = 0;
kv = 0;
while  k < m
    if k
        q = A(v_next);
    else
        q = v_start;
    end
    
    % Orthogonalize q against previous vr's
    qc = orthagainst([qr;qi],Vr(:,1:k));
    k = k + 1;
    Vr(:,k) = qc;
    v_next = complex(qc(1:n),qc(n+1:end));
    
    V2(:,2*k-1:2*k) = [qc(1:n) qc(n+1:end)];
    
%     % Orthogonalize [Re(q) Im(q)] against previous v's
%     q = orthagainst([qr qi],V(:,1:kv));
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