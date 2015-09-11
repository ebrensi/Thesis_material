function [Ar Rr0] = makeAr(C,G,b,s0)
% [Ar Rr0] = makeAr(C,G,b,s0)
% returns a handle to the function that implements the matrix-vector 
%  multiplication Ar*v (equivalent real) that we use in the Arnoldi process. 

N = length(b);

C_s0G = C-s0*G;
if issparse(C)
    [L U P Q] = lu(C_s0G);
    Ar = @multAr;
else
    [L,U,P] = lu(C_s0G);
    Ar = @multAr_dense;
end

if nargout == 2
    r0 = -C_s0G\b;
    Rr0 = [real(r0); imag(r0)];
end

% -----------------------------------------------
    function qr = multAr(vr)
        v = complex(vr(1:N,:),vr(N+1:end,:));
        y = P*(G*v);
        q = Q * (U\(L\y));
        qr = [real(q); imag(q)];
    end

    function qr = multAr_dense(vr)
        v = complex(vr(1:N,:),vr(N+1:end,:));
        y = P*(G*v);
        q = U\(L\y);
        qr = [real(q); imag(q)];
    end
end
