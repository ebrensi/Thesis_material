function tfunc = transfer_function(V,A,E,C,B,s)
% tfunc = transfer_function(V,A,E,C,B,s)
%
%  Evaluates the transfer function
%    H(s) = Cn^T (An-sEn) \ Bn
%  at point(s) s.
%
%  where (An,En,Cn,Bn) is the realization (A,E,C,B) orthogonally projected onto V.
%  if V = [] then we assume V = I.


m = size(B,2); p = size(C,2);  % p, m  are input, output dimensions
siso = m == 1 & p == 1;  % is this a SISO model?

Ls = length(s(:));


if isempty(V)
	An = A; En = E; Cn = C; Bn = B;
else
	An = V'*A*V;
	En = V'*E*V;
	Bn = V'*B;
	Cn = V'*C;
end

if issparse(An) && issparse(En)
	if siso
		tfunc = zeros(1,Ls);
		for j = 1:Ls
			[L U P Q] = lu((s(j)*En-An));
			tfunc(j) = Cn'*(Q*(U\(L\(P*Bn))));
		end
	else
		tfunc = zeros(m,p,Ls);
		for j = 1:Ls
			[L U P Q] = lu((s(j)*En-An));
			tfunc(:,:,j) = Cn'*(Q*(U\(L\(P*Bn))));
		end
	end
	
else
	%  This is for non-sparse realization 
	[TAn,TEn,Q,Z] = qz(An,En);
	Cn = Cn'*Z;
	Bn = Q*Bn;
% 	optU.UT = true;
	if siso
		tfunc = zeros(1,Ls);
		for j = 1:Ls
            tfunc(j) = Cn *( (s(j)*TEn-TAn) \ Bn ) ;
		end
	else
		tfunc = zeros(m,p,Ls);
		for j = 1:Ls
			tfunc(:,:,j) = Cn*((s(j)*TEn-TAn)\Bn);
		end
	end
end

if siso
	tfunc = reshape(tfunc,size(s));
end
