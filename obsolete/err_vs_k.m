function err = err_vs_k(data,kmax)
% outputs err for k = 1:kmax

[C G c b] = inputdata(data);
frq = logspace(8,10,200);  
s = 2*pi*sqrt(-1)*frq;

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);

s0 = pi * 1e10;
[A r0] = makeA(C,G,b,s0);

err = step_through_arnoldi(A,r0,kmax,@tfunc_err,2);
[min_err ix] = min(err);

semilogy(1:kmax,err,'.');%,ix,min_err,'r+');
legend('hess','PRIMA');
title(sprintf('%s: relerr vs k',data));

%% ---------------------------------------------------------
    function err = tfunc_err(Vk,Hk)
        tfunc = abs(tf_hess(Vk,Hk,r0,c,s,s0));
        err_hess = norm(tfunc_real(:) - tfunc(:),inf)/ntf;
        
        tfunc = abs(tf_proj(Vk,C,G,c,b,s));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;
        
        err = [err_hess err_proj];
    end
end
