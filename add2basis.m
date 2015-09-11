function [Y T kept] = add2basis(Y, T, y,tol)
% [Y T kept] = add2basis(Y, T, y, tol)
%
% Assuming (Y,T) is a QR decomposition, we add y to it with inexact deflation.
ell = size(Y,2);
k = size(y,2);
kept = false(1,k);

for i = 1:k
	t = zeros(ell+1,1);
	for j = 1:ell
		t(j) = Y(:,j)' * y(:,i);
		y(:,i) = y(:,i) - t(j)*Y(:,j);  % subtract off directions Y.
	end
	t(ell+1) = norm(y(:,i));
	if t(ell+1) > tol
		kept(i) = true;
		Y = [Y y(:,i)/t(ell+1)];
		T = [[T; zeros(1,ell)] t];   
		ell = ell + 1;
	else
		kept(i) = false;
% 		fprintf('Ritz-vec deflated: %g < %g\n',t(ell+1),tol)
	end
end
