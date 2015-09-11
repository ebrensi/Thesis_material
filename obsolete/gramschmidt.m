function [B count] = gramschmidt(A,tol)
% [B count] = gramschmidt(A,tol)
%
% B is the orthnormalization of A
% count(i) is the number of columns of B that span columns 1:i of A. 
%  ex. B(:,1:count(4)) spans the same space as A(:,1:4)

% QR decomposition
[Q,R]=qr(A,0);
% Determine effective rank
if ~exist('tol','var') || isempty(tol)
    tol = eps*norm(A,'fro');
end
LI = abs(diag(R)) > tol;
count = cumsum(LI);
B = A(:,LI);
