function [mu wt] = p1(A,E,c,b)
[AA,BB,Q,Z,V,W] = qz(full(A),full(E),'complex');
aa = diag(AA);
bb = diag(BB);
mu = aa ./ bb;


f = V'*c;
g = W'*b;

x = conj(f) .* g;

e = diag(W'*E*V);

wt = abs(x) ./ abs(e) ./ abs(real(mu));  % weight as implied by Joost Rommes

[wt idx] = sort(wt);
mu = mu(idx);

end
