function [AE R0] = make_SI_op(A,E,b,s0)
% [AE R0] = make_SI_op(A,E,b,s0)
%
% returns a handle to the shift-invert operator AE(s0) and start "residual" R0
%   that we use in the rational-Arnoldi process. 

if isinf(s0)
    % multiplier and start vector for infinite expansion point 
    
else
    A_s0E = A-s0*E;
    if issparse(A_s0E)
        [L U P Q] = lu(A_s0E);
        AE = @multAE;
    else
        [L,U,P] = lu(A_s0E);
        AE = @multAE_dense;
    end
    
    if nargout > 1
        R0 = -A_s0E\b;
    end
end
% -----------------------------------------------
    function q = multAE(v)
        y = P*(E*v);
        q = Q * (U\(L\y));
    end

    function q = multAE_dense(v)
        y = P*(E*v);
        q = U\(L\y);
    end
end
