function [mu wt] = p2(C,G,c,b)
s0 = pi * 1e10;

%% ****** compute eigenvalues of A(s0) = (C-s0*G)^-1 * G  ******
s0G_C = s0*G-C;
r0 = s0G_C\b;
% ** compute eigs of A **
A = -s0G_C\G;
[Z L] = eig(full(A));
L = diag(L);
mu = s0 + 1./L;

%% compute transfer function terms
f = Z'*c;
g = Z\r0;

fg = conj(f) .* g;
wt = abs(fg) .* abs(s0-mu) ./ abs(real(mu));

[wt idx] = sort(wt);
mu = mu(idx);
end
