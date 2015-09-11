function pr = rmpassive(data,kvals,s0,meth)
% pr = rmpassive(data,kvals,s0,meth)

[C G c b] = inputdata(data);
[A r0] = makeA(C,G,b,s0);

K = kvals(end);
bk_hess = [norm(r0); zeros(K-1,1)];

if all(b == c)
    do_alt = true;
else
    do_alt = false;
end
    
if strcmpi(meth,'prima')
    pr = step_through_arnoldi(A,r0,kvals,@PRIMA_pcheck,1);
else
    pr = step_through_arnoldi(A,r0,kvals,@hessenberg_pcheck,1);
end

% ------------------------------------------------------------------
    function tf = PRIMA_pcheck(Vk,Hk)
        % make the PRIMA reduced model
        %  i.e. H(s) = ck' * (sGk - Ck)^-1 * bk
        if ~isreal(Vk)
            Vk = [real(Vk) imag(Vk)];
        end
        Ck = Vk'*C*Vk;
        Gk = Vk'*G*Vk;
        bk = Vk'*b;
        ck = Vk'*c;

        % determine whether it is positive real
        tf = isposreal(Ck,Gk,ck,bk,s0);
        if do_alt
            tf2 = isposreal2(Ck,Gk);
            fprintf('%d: [tf tf2]: %d %d\n',size(Vk,2),tf,tf2);
        else
            fprintf('%d: %d\n',size(Vk,2),tf);
        end
    end

    function tf = hessenberg_pcheck(Vk,Hk)
        % make the reduced model via Arnoldi generated Hessenberg matrix
        %  i.e. H(s) = ck' * (I - (s-s0)Hk)^-1 * bk
        %            = -ck' * (sHk - (I+s0Hk) )^-1 * bk
        k = length(Hk);
        bk = bk_hess(1:k);
        ck = Vk'*c;
        Ik = eye(k);
        Is0Hk = Ik+s0*Hk;
        tf = isposreal(Is0Hk, Hk,-ck,bk,s0);
        fprintf('%d: %d\n',size(Vk,2),tf);
    end
end
