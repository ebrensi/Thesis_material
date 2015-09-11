function [mu W rr ritz_vals tf_terms] = approx_poles_hess(V,H,vnext,c,rho,s0)
% [mu Z rr ritz_vals tf_terms] = approx_poles_hess(V,H,vnext,c,rho,s0)
% 
%  Compute the pole/eigenvalue decomposition of a ROM from its Arnoldi
%  decomposition, via eigenvalues of the Hessenberg matrix.


s = getS();
SISO = (size(c,2) * size(rho,2) == 1);

[W D] = eig(H);
ritz_vals = diag(D);

rr = abs(W(end,:)*norm(vnext))'./abs(ritz_vals); % relative residual 2-norm
mu = s0 + 1./ritz_vals; % poles

if any(isinf(mu))
   warning('infinite pole') 
end

if nargout > 4
    % compute transfer function terms
    f = c'*V*W;
    g = W \ rho;
    if SISO
        x = f' .* g;
    else
        Lmu = length(mu);
        x = zeros(Lmu,1);
        for k = 1:Lmu
            Xk = f(k,:)'*g(k,:);
            x(k) = norm(Xk,inf);
        end
    end
    
    tf_terms = bsxfun(@rdivide, x, 1-ritz_vals*(s-s0));
end

end
