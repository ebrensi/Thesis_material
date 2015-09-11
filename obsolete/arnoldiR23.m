function [basis_splitC basis_splitR VsplitC VsplitR ] = arnoldiR23(A,r,m)
% V = arnoldiR(A,v_start,m)
%
%  This is the real-Arnoldi process for computing a basis for the k-th 
%  Split-Krylov subspace K_k(A,r)*.

if isa(A,'numeric')
    A = @(x) A*x;
end

n = length(r);
v_start = r;
rclass = class(r);

VC = zeros(n,m,rclass);
VR = zeros(2*n,m,rclass);

VsplitC = zeros(n,2*m,rclass);
VsplitR = zeros(n,2*m,rclass);

k = 0;
while  k < m
    if k
        qC = A(vC_next);
        qR = A(vR_next);
    else
        qC = v_start;
        qR = v_start;
    end

    % Orthogonalize q against previous vc's
    vC_next = orthagainst(qC,VC(:,1:k));
    qc = orthagainst([real(qR);imag(qR)],VR(:,1:k));
    
    k = k + 1;
    VC(:,k) = vC_next;
    VsplitC(:,2*k-1:2*k) = [real(vC_next) imag(vC_next)];
     
    VR(:,k) = qc;
    qc1 = qc(1:n);  qc2 = qc(n+1:end);
    vR_next = complex(qc1,qc2);
    VsplitR(:,2*k-1:2*k) = [qc1 qc2];
end

% orth is a matlab built-in function that uses the orthogonal 
%  basis produced via SVD 
%  orthagainst uses grahm-schmidt 
[basis_splitC] = orth(VsplitC);
[basis_splitR]= orth(VsplitR);