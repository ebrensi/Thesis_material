function Q = my_gs(A,tol)
%
[N M] = size(A);
%
Q = zeros(N,0);
%
k = 1;
%
for j = 1:M,
  q = A(:,j);
  for i = 1:k-1,
    alpha = (Q(:,i))' * q;
    q = q - alpha * Q(:,i);
  end
  normq = norm(q);
  if normq >= tol,
    Q(:,k) = q / normq;
    k = k + 1;
  end
end  
