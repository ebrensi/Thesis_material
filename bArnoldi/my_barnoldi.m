function [V,H] = my_barnoldi(MultA,R,nmax)
%
%
%   Block Arnoldi procedure 
% 
%   This version uses modified Gram-Schmidt and
%   no deflation
%   (as in the PRIMA paper)
%
%   m = number of block steps
% 
%   B = starting block
%
%   Matrix multiplications are done via the function MultA
%
%initialization and normalization

[N,m] = size(full(R));

[Q,RR]=qr(full(R),0);

V(:,1:m) = Q;

%
%   Arnoldi loop (Gram-Schmidt orthogonalization process) 
%   A*V = V*H + ( 0,....,0, Z)  

j=1; 
while j <= nmax, 
%
   jl = (j-1)*m+1;
   ju = j*m;
   Z  = feval(MultA,V(:,jl:ju));
%   
   for i = 1:j
      il = (i-1)*m+1;
      iu = i*m;      
      Hij = V(:,il:iu)'*Z;
      H(il:iu,jl:ju) = Hij;
      Z = Z - V(:,il:iu)*Hij;
   end
%    
   [Q,RR]=qr(full(Z),0);
   V(:,jl+m:ju+m) = Q;
   H(jl+m:ju+m,jl:ju) = RR;
%
   j = j+1;
end
%
