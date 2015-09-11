function [Vr W] = make_basis_real(V)
% [Vr W] = make_basis_real(V)

% Constructs a real basis Vr that spans the same space as V, in the same order.

if ~isreal(V)
	N = size(V,1);
	n = size(V,2);
	W = zeros(N,2*n);
	
	k = 1;
	for j = 1:n
		vj = V(:,j);
		W(:,k) = real(vj);
		W(:,k+1) = imag(vj);
		k = k+2;
	end
	
% 	Orthogonalize
	[Vr R] = qr(W,0);
else 
	Vr = V;
	W = V;
end
