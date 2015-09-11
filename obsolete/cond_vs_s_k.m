function opt = cond_vs_s_k(data,kmax)
% outputs cond(Ck-sGk) and (Ik-(s-s0)Hk) for varying s and k=1:kmax

[C G c b] = inputdata(data);
frq = logspace(8,10,200);  
s = 2*pi*sqrt(-1)*frq;
Ls = length(s);

s0 = pi * 1e10;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

opt = step_through_arnoldi(@A,r0,kmax,@ccg,Ls);

%% ---------------------------------------------------------
    function cc = ccg(Vk,Hk,k)
%         Ck = Vk'*C*Vk;
%         Gk = Vk'*G*Vk;
        cc = zeros(1,Ls);
        I = eye(k);
        for j = 1:Ls
%            cc(j) = cond(s(j)*Gk-Ck);
           cc(j) = cond(I-(s(j)-s0)*Hk);
        end
        k
    end

    function q = A(v)
        q = C_s0G\(G*v);
    end

end
