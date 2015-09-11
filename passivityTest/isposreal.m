function posreal = isposreal(A,E,c,b,s0)
% posreal = isposreal(A,E,c,b,s0)
%
% this function implements a test for whether an RCL circuit model is
% passive, devised by Z. Bai & R. W. Freund, 12/2000 in
% "Eigenvalue-Based Characterization and test for Positive Realness
%  of Scalar Transfer Functions"
%  transfer fuction given by H(s) = c' * (sE-A)^-1 * b

verbose = true;
rel_diff_tol = 0.2;
zero_tol = 0.01;

n = size(A,1);
I = speye(n);

[gamma S T] = reduce2unitvecs(A,E,c,b,s0);

if isnan(gamma)  % i.e. if E=I and gamma=0
    posreal = false;
    return;
end

%%  Step 2. of algorithm 2
S = full(S);  T = full(T);
% poles (mu) of the transfer function are the eigenvalues of the pencil (S,T)
mu_ST = eig(S,T);
[mu multiplicity cnj1] = process_eigs(mu_ST);
% dbug1 = [mu multiplicity cnj1];

if isinf(mu(end))
    mu(end) = [];                         % discard infinite poles
    multiplicity(end) = [];
end

pure_imag_mu = real(mu)==0 & imag(mu)~=0;
pos_real_ind = real(mu)>0;

if any(pos_real_ind) || any(multiplicity(pure_imag_mu)>1)
    % (or if any pure_imag pole has non-positive residue)
    %  this still needs to be implemented
    S_hat = S(2:end,2:end);  T_hat = T(2:end,2:end);
    [mu_hat mul_hat cnj2] = process_eigs(eig(S_hat,T_hat));
%     dbug2 = [mu_hat mul_hat cnj2];
    if isinf(mu_hat(end))
        mu_hat(end) = [];                         % discard infinite poles
        mul_hat(end) = [];
    end
    pos_real_hat_ind = real(mu_hat)>0;

    % (S,T) and (S_hat,T_hat) should have the same eigs w/ pos real part,
    % to some tolerance
    num_pos_real = nnz(multiplicity(pos_real_ind));
    num_pos_real_hat = nnz(mul_hat(pos_real_hat_ind));

    pos_real_mu = expand_elements(mu(pos_real_ind),multiplicity(pos_real_ind));
    pos_real_mu_hat = expand_elements(mu_hat(pos_real_hat_ind),mul_hat(pos_real_hat_ind));

    eigsame =  num_pos_real == num_pos_real_hat;
    if eigsame
        % make sure they're the same ones (note that at this point they are sorted)
        pos_real_mu = expand_elements(mu(pos_real_ind),multiplicity(pos_real_ind));
        pos_real_mu_hat = expand_elements(mu_hat(pos_real_hat_ind),mul_hat(pos_real_hat_ind));
        chk = [pos_real_mu pos_real_mu_hat];
        rel_diff_between_posreals = norm(pos_real_mu-pos_real_mu_hat)/norm(pos_real_mu);
        
        eigsame = rel_diff_between_posreals < rel_diff_tol;
    end
    if ~eigsame
        if verbose
            fprintf('pos real eigs of (S_hat,T_hat) differ from those of (S,T)\n');
        end
        posreal = false;
        return;
    end
    
end

%% step 3.
e1 = eye(n,1);
M = S^2;   M(1:end,1) = S(1:end,1); % i.e.  M = S^2 *(I - e1e1T) + S*e1e1T;
N = T^2;   N(1:end,1) = 0;          % i.e.  N = T^2 *(I - e1e1T);
tau = norm(M,inf)/norm(N,inf);  % Scaling factor for better numerical accuracy.
[AA BB] = qz(M,-tau*N,'complex');
aa = diag(AA);
bb = diag(BB);

aa_small = magvals(aa,1);
bb_small = magvals(bb,1);

% if the pencil M+sN is singular then H(s) is lossless positive real
if any(aa_small & bb_small)
    posreal = true;
    if verbose 
        fprintf('H(s) is lossless positive.\n');
    end
    return;
end

% here we may assume that M+sN is not singular
MN_eigs = aa./bb;
[lambda multiplicity] = process_eigs(MN_eigs);
lambda(isinf(lambda)) = [];
posreal = real(lambda)>0;
if any(isodd(multiplicity(posreal)));
    % if any positive real lambda has odd multiplicity
    posreal = false;
    if verbose
        fprintf('(M,N) has a pos real eigenvalue with odd multiplicty\n')
    end
    return;
end

%% setp 6.
% we pick lambda0 so that det(M + lambda0*N) ~= 0
lambda0 = ceil(max(abs(lambda))*2);  
% lambda0 = tau;

%matlab recommends against using det to test for singularity
d = gamma*det(M+lambda0*N);  
if d < -eps
    posreal = false;
    if verbose
        fprintf('step 6 test fails\n')
    end
    return;
end

%% tfunc meets all criteria
posreal = true;


%% --------------------------------------------------------------------

% this is algorithm 1 from the paper. returns T,S such that ST = TS and
%  H(s) = gamma*e_1'* (S - sT) \ e_1 = c' * (sE - A)\b
% i.e. Transfer function formulated as multiplication with unit vectors.
    function [gamma S T] = reduce2unitvecs(A,E,c,b,s0)
        if norm(E-I,inf) < eps  % i.e. if E = I
            gamma = -c'*b;
            if abs(gamma) < eps
                gamma = nan;
                return;
            else
                v = -b;   w = c;
                S = lemma1(v,w,A);
                T = I;
            end
        else                        %i.e. if E ~= I
            As0E = A-s0*E;
            if issparse(As0E)
                [L U P Q] = lu(As0E);
            else
                [L U P] = lu(As0E);
                Q = speye(n);
            end

            gamma = -c' * Q * (U\(L\(P*b)));  % gamma = c'*(s0*E-A)\b;
            G = L\((P*E*Q)/U);
            v = -L\(P*b);    w = U' \ (Q'*c);

            T = lemma1(v,w,G);
            S = I + s0*T;
        end

        % ----------------------------------------
        function C = lemma1(v,w,M)
            % compute  result = Q*M*inv(Q) where Q is a theoretical matrix
            % (not computed explicitly) that satifies (7) from the paper.
            v = v(:); w = w(:);
            e1 = eye(n,1);
            [q1,beta1,gamma_v] = gallery('house',v);

            Q1w = w - beta1*(q1'*w)*q1;
            gamma_w = Q1w(1);
            w_tilde = Q1w(2:end);

            [q2,beta2,eta] = gallery('house',w_tilde);

            B = HHproj(q1,beta1,M);
            C = HHproj([0; q2],beta2,B);
            C1 = C(1:2,1:2);
            C2 = C(1:2,3:end);
            C3 = C(3:end,1:2);
            R = [1 eta/gamma_w; 0 1];
            Rinv = eye(2); Rinv(1,2) = -R(1,2);

            % Explicitly compute Q and Q2 just to show that this works
            %             Q1 = I - beta1*q1*q1';
            %             Q2 = I - beta2*[0;q2]*[0;q2]';
            %             RR = eye(n); RR(1:2,1:2) = R;
            %             Q = RR*Q2*Q1;
            %             Calt = Q*M*inv(Q);

            C(1:2,1:2) = R*C1*Rinv;
            C(1:2,3:end) = R*C2;
            C(3:end,1:2) = C3*Rinv;
            % ---------------------------------
            function PAP = HHproj(v,b,A)
                % compute P*A*P, where P=I-b*v*v'
                AP = A - b*(A*v)*v';
                PAP = AP - b*v*(v'*AP);
            end
        end
    end


    function out = expand_elements(elements,multiplicity)
        out = [];
        Lin = length(elements);
        for i = 1:Lin
            e = elements(i);
            out = [out; e(ones(multiplicity(i),1))];
        end
    end

    function result = isodd(mat)
        % returns a logical array indicating which elements of mat are odd;
        result = ~mod(mat,2);
    end

end
