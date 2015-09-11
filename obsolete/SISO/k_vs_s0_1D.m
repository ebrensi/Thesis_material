function [min_K s0_opt] = k_vs_s0_1D(data,nkpts,ns0pts)
% [min_K s0_opt] = k_vs_s0_1D(data,nkpts,ns0pts)
%
% determines number of steps k<kmax required to achieve a tfunc error of err_tol
%  for several values of s0

kmax = 200;
kvals = fix(linspace(2,kmax,nkpts));
err_tol = 0.01;

% rng = linspace(1e-6,3,ns0pts);
rng = linspace(1e-6,2.5,ns0pts);
S0 = rng*sqrt(-1)*pi*1e10;
% S0 = rng*pi * 1e10;

[C G c b] = inputdata(data);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);
K = zeros(ns0pts,1);

tic;
for j = 1:ns0pts
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    err = step_through_arnoldi(A,r0,kvals,@tfunc_err,1);
    k_proj = find(err);
    if isempty(k_proj) 
        k_proj = kmax;
    else
        k_proj = kvals(k_proj);
    end
    K(j) = k_proj;
    et = toc;
    rt = et/j;
    eta = rt*(ns0pts-j);
    fprintf('eta: %d:%02d\r',fix(eta/60),fix(rem(eta,60)));
end
fprintf('\n');

[min_K idx] = min(K);
s0_opt = S0(idx);
plot(rng,K,'-',rng(idx),min_K,'r+');
title(sprintf('%s:\t tol = %2.1e',data,err_tol));
legend(sprintf('proj: (%4.3g, %d)',rng(idx),min_K));

saveas(gcf,sprintf('%s_k_vs_s0_1D',data));
% -----------------------------------------------------
    function errout = tfunc_err(Vk,Hk)
        tfunc = abs(tf_proj(Vk,C,G,c,b));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;

        if err_proj < err_tol 
            errout = nan;
        else
            errout = err_proj < err_tol;
        end
    end

end  % end main function
