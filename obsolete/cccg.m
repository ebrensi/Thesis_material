function [cC cG] = cccg(data,kmax)
% outputs cond(Ck) and cond(Gk) for k = 1:kmax

[C G c b] = inputdata(data);
frq = logspace(8,10,200);  
s = 2*pi*sqrt(-1)*frq;

s0 = pi * 1e10;
C_s0G = C-s0*G;
r0 = -C_s0G\b;

opt = step_through_arnoldi(@A,r0,kmax,@ccg,2);
cC = opt(:,1); 
cG = opt(:,2);

%% ---------------------------------------------------------
    function opt = ccg(Vk,Hk,k)
        Ck = Vk'*C*Vk;
        Gk = Vk'*G*Vk;
        opt = [cond(Ck) cond(Gk)];
    end

    function q = A(v)
        q = C_s0G\(G*v);
    end

end
