%  This is a code snippit to add schur decomp functionality to the deflate
%  function.


[X U] = schur(H,'complex');   % compute full Schur decomposition of H
[L idx] = sort(ordeig(U));    % sort eigs by magnitude
[X U] = ordschur(X,U,cvrgd(idx));  % re-order schur decomp
X = X(:,1:n_cvrgd);
U = U(1:n_cvrgd,1:n_cvrgd);
L = diag(U);
L_next = L ./ (1-(s0_next-s0)*L);
U_next = U;
U_next(1:(n_cvrgd+1):end) = L_next;

%             chk1 = [sort(lambda(cvrgd))  sort(L)]
%             chk = [sort(lambda_next(cvrgd)) sort(L_next)];

uT = X(end,:);  % the last row of Q
Up = [U_next; eta*uT];
Y = [Yold V]*X;
