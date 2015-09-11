function tft = tfterms(lambda,Z,c,r0,s0)
% [tft wt] = tfterms(lambda,Z,c,r0)
%
%  returns the terms in the pole-residue formulation of the transfer
%  function given by H(s) = c'*(A-sE)\b = c'*(I-(s-s0)H)\r0

SISO = (size(c,2) * size(r0,2) == 1);

f = c'*Z;
g = Z \ r0;
if SISO
    x = f' .* g;
else
    Lmu = length(lambda);
    x = zeros(Lmu,1);
    for k = 1:Lmu
        Xk = f(k,:)'*g(k,:);
        x(k) = norm(Xk,inf);
    end
end

s = getS();
tft = bsxfun(@rdivide, x, 1-lambda*(s-s0));

