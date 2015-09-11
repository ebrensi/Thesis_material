function [mu Z rr tf_terms] = approx_poles_proj(V,A,E,C,B,s0,cconj,meth)
% [mu Z rr tf_terms] = approx_poles_proj(V,A,E,C,B,s0)
%
%  Given H(s) = C' * (sE-A) * B,
%  Computes the pole/eigenvalue decomposition of H(s) from its Arnoldi
%  decomposition, via eigenvalues of the projected pencil (An,En).

rr_tol = eps;
sc = size(C,2);
sb = size(B,2);
SISO = (sc * sb == 1);

s = getS;

if isempty(V)
    An = A; En = E; Bn = B; Cn = C;
else
    An = V'*A*V;
    En = V'*E*V;
    Bn = V'*B;
    % we expect B == C in most cases
    if any(B ~= C)
        Cn = V'*C;
    else
        Cn = Bn;
    end
end

% for better numerical results, we introduce a scaling to 
%   ensure that the norms of An and En are of the same order.
%  we do this via floating point exponent manipulation 
% [f1 eA] = log2(norm(An,inf));
% [f2 eE] = log2(norm(En,inf));
% ediff = eA-eE;
% t1En = pow2(En,ediff);  % shift exponents of elements in En to match those of An

% tau = norm(An,inf)/norm(En,inf);
% t2En = tau*En;

if nargout > 3
    if ~exist('meth','var') || strcmp(meth,'e')
        [mu W tf_terms] = urm_e(An,En,Cn,Bn);
    elseif strcmp(meth,'p')
        [mu W tf_terms] = urm_p(An,En,Cn,Bn);
    elseif strcmp(meth,'qz')
        [mu W tf_terms] = pole_decomp(An,En,Cn,Bn);
    end
else
    if ~exist('meth','var') || strcmp(meth,'e')
        [mu W] = urm_e(An,En,Cn,Bn);
    elseif strcmp(meth,'p')
        [mu W] = urm_p(An,En,Cn,Bn);
    elseif strcmp(meth,'qz')
        [mu W] = pole_decomp(An,En,Cn,Bn);
    end
end

% Shift elements of mu to their unscaled values
% rmu = pow2(real(mu_t1),ediff);
% imu = pow2(imag(mu_t1),ediff);
% mu1 = complex(rmu,imu);  
% 
% % mu = tau*mu_t2;
% mu_pow2 = mu1; 
% mu_tau = mu;
% mu_straight = mu_2;


Z = V*W;
M = diag(mu);
AZ = A*Z;

% compute relative residuals of the ritz values
rr = max(abs(AZ - E*Z*M))' ./ max(abs(AZ))';

if any(isinf(mu)) 
   warning('infinite pole') 
end

if exist('cconj','var') && ~isempty(cconj)
    % combine tf_terms for conjugate pairs and hold on to only one half
    ess = imag(mu) >= 0;
    mu_ess = mu(ess);
    cplx = imag(mu_ess) > 0;
    tf_terms_ess = tf_terms(ess,:);
    tf_terms_ess(cplx,:) = tf_terms_ess(cplx,:) + tf_terms(~ess,:);
    tf_terms = tf_terms_ess;
    rr = rr(ess);
    mu = mu_ess;
end


%% ***************************************************************
    function [mu Z tf_terms] = urm_p(A,E,c,b)
        Ls = length(s);
        %% ** compute poles as eigenvalues mu of the pencil (sE-A)
        % eig(A,E) is the correct formulation but it returns inf eigenvalues
        %  which cause trouble with matlab. So we use eig(E,A) and get
        %  eigenvalues 1/mu.
        [Z invM] = eig(full(E),full(A));
        inv_mu = diag(invM);
        mu = 1 ./ inv_mu;
        Lp = length(mu);

        %% *** compute transfer function terms
        if nargout > 2
            f = Z'*c;
            g = ((s0*E-A)*Z)\b;
            if SISO                
                x = conj(f) .* g;
            else 
                Lmu = length(mu);
                x = zeros(Lmu,1);
                for k = 1:Lmu
                    Xk = f(k,:)'*g(k,:);
                    x(k) = norm(Xk,inf);
                end
            end
            finite_mu = abs(inv_mu) ~= 0;
            zero_mu = isinf(inv_mu);
            fnz_mu = finite_mu & ~zero_mu;
            s0_mu = s0*inv_mu(fnz_mu) - 1;
            s_mu = bsxfun(@times,s,inv_mu(fnz_mu)) - 1;
            tf_terms = zeros(Lp,Ls);
            tf_terms(~finite_mu,:) = repmat(x(~finite_mu),1,Ls);
            tf_terms(zero_mu,:) = s0*bsxfun(@rdivide,x(zero_mu),s);
            tf_terms(fnz_mu,:) = bsxfun(@rdivide,x(fnz_mu).*s0_mu,s_mu);
        end
    end

%%
    function [mu Z tf_terms] = urm_e(A,E,c,b)
%         [f1 eA] = log2(norm(A,inf));
%         [f2 eE] = log2(norm(E,inf));
%         [f3 es0] = log2(norm(s0));
%         ediff = eA-eE-es0;
%         ediff = 0;
%         tE = pow2(E,ediff);  % shift exponents of elements in En to match those of An
        tE = E;
        
        %% ****** compute eigenvalues of (A-s0*E)^-1 * E  ******
        s0tE_A = s0*tE-A;
        r0 = s0tE_A\b;
        % ** compute eigs of A **
        M = -s0tE_A\tE;
        [Z L] = eig(full(M));
        L = diag(L);
        mu = s0 + 1./L;

        if nargout > 2
            %% compute transfer function terms
            f = Z'*c;
            g = Z\r0;
            if SISO                
                x = conj(f) .* g;
            else 
                Lmu = length(mu);
                x = zeros(Lmu,1);
                for k = 1:Lmu
                    Xk = f(k,:)'*g(k,:);
                    x(k) = norm(Xk,inf);
                end
            end
            
            tf_terms = bsxfun(@rdivide, x, 1-L*(s-s0));
        end
        
        % Shift elements of mu to their unscaled values
%         rmu = pow2(real(mu),ediff);
%         imu = pow2(imag(mu),ediff);
%         mu = complex(rmu,imu);
    end

%%
    function [mu V tf_terms] = pole_decomp(A,E,c,b)
        % for better numerical results, we introduce a scaling to
        %   ensure that the norms of A and E are of the same order.
        %  we do this via floating point exponent manipulation
%         [f1 eA] = log2(norm(A,inf));
%         [f2 eE] = log2(norm(E,inf));
%         ediff = eA-eE;
        ediff = 0;
        tE = E;
%         tE = pow2(E,ediff);  % shift exponents of elements in En to match those of An

        [AA,BB,Q,Z,V,W] = qz(full(A),full(tE),'complex');
        aa = diag(AA);
        bb = diag(BB);
        mu = aa ./ bb;
        
        if nargout > 2
            f = V'*c;
            g = W'*b;
            if SISO                
                x = conj(f) .* g;
            else 
                Lmu = length(mu);
                x = zeros(Lmu,1);
                for k = 1:Lmu
                    Xk = f(k,:)'*g(k,:);
                    x(k) = norm(Xk,inf);
                end
            end
            
            e = diag(W'*E*V);
            a = diag(W'*A*V);
            se_a = bsxfun(@minus,e*s,a);
            tf_terms = bsxfun(@rdivide, x, se_a);
        end
       
        % Shift elements of mu to their unscaled values
%         rmu = pow2(real(mu),ediff);
%         imu = pow2(imag(mu),ediff);
%         mu = complex(rmu,imu);
    end
end
