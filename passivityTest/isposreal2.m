function posreal = isposreal2(A,E)
% posreal = isposreal2(A,E)
%
% This is the test for positive realness of the transfer function
% H(s) = b' * (sE-A)^1 * b (note both left and right side b),
% implied by theorem 4.5 of
% "Model reduction methods based on Krylov subspaces", R. W. Freund
% Acta Numerica (2003)

% we need to check if A+A' is negative semi definite and E is positive semi-definite.
% and the pencil sE-A is regular (non-singular)

posreal = false;
if ispsd(-(A+A'))
    if ispsd(E);
        if ~issingular(A,E)
            posreal = true;
        else
            v = '(A,E) is singular';
        end
    else
        v= 'E is not positive semidefinite';
    end
else
    v = sprintf('A\''+ A is not negative semidefinite');
end

if ~posreal
    fprintf('%s\n',v);
end

%% --------------------------------------------------------
    function psd = ispsd(M)
        dM = diag(M);
        maxM = max(max(abs(M)));
        if any(dM < -eps) || maxM ~= max(dM)
            psd = false;
            return;
        end

        zd_idx = find(dM==0);
        if any(any(M(zd_idx,:))) || any(any(M(:,zd_idx)))
            psd = false;
            return;
        end

        % now we may assume necessary conditions are met
        eigM = eig(full(M));
        lambda = process_eigs(eigM);
        if any(lambda < 0)
            psd = false;
        else
            psd = true;
        end

    end

end  % main function 
