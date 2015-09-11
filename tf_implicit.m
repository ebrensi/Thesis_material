function tfunc = tf_implicit(H,rho,c,s0,s)
% tfunc = tf_implicit(H,rho,c,s0,s)
%
% obtaining the transfer function via hessenberg Matrix H
% H(s) = (V'*c)' * (I - (s-s0)H)^-1 * (V'*r0) 
%   note that we expect rho = V'*r0 to be input

nout = size(c,2); nin = size(rho,2);
siso = nout == 1 & nin == 1;  % is this a SISO model?

n = length(H);
I = eye(n);
Ls = length(s(:));

if siso
		tfunc = zeros(1,Ls);
		for j = 1:Ls
			tfunc(j) = c' * ( (I - (s(j)-s0)*H)\rho ) ;
		end
	else
		tfunc = zeros(nout,nin,Ls);
		for j = 1:Ls
			tfunc(:,:,j) =  c' * ( (I - (s(j)-s0)*H)\rho ) ;
		end
end

if siso
	tfunc = reshape(tfunc,size(s));
end
