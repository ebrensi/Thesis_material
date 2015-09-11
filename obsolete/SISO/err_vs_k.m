function err = err_vs_k(data,kmax,npts,s0)
% outputs err for k = 1:kmax

[C G c b] = inputdata(data,true);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);

% s0 = pi * 1e10;
[A r0] = makeA(C,G,b,s0);

kvals = fix(linspace(2,kmax,npts));
err = step_through_arnoldi(A,r0,kvals,@tfunc_err,2);
[min_err ix] = min(err);

semilogy(kvals,err,'.',kvals(ix),min_err,'r+');
legend('hess','PRIMA');
title(sprintf('%s: relerr vs k  \t   at s0 = %3.2g + %3.2gi',data,real(s0),imag(s0)));

%% ---------------------------------------------------------
    function err = tfunc_err(Vk,Hk)
        tfunc = abs(tf_hess(Vk,Hk,r0,c,s0));
        err_hess = norm(tfunc_real(:) - tfunc(:),inf)/ntf;
        
        tfunc = abs(tf_proj(Vk,C,G,c,b));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;
        
        err = [err_hess err_proj];
    end
end
