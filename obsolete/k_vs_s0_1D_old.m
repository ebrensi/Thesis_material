function [min_K s0_opt] = k_vs_s0_1D(data,nkpts,ns0pts)
% [min_K s0_opt] = k_vs_s0_1D(data,nkpts,ns0pts)
%
% determines number of steps k<kmax required to achieve a tfunc error of err_tol
%  for several values of s0

kmax = 200;
kvals = fix(linspace(2,kmax,nkpts));
err_tol = 0.01;

rl_rng = linspace(1e-6,2,ns0pts);
S0 = rl_rng *pi * 1e10;

[C G c b] = inputdata(data);
frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);
K = zeros(ns0pts,2);

tic;
for j = 1:ns0pts
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    err = step_through_arnoldi(A,r0,kvals,@tfunc_err,1);
    k_hess = find(err(:,1),1);
    if isempty(k_hess) 
        k_hess = kmax;
    else
        k_hess = kvals(k_hess);
    end
    k_proj = find(err(:,2),1);
    if isempty(k_proj) 
        k_proj = kmax;
    else
        k_proj = kvals(k_proj);
    end
    K(j,:) = [k_hess k_proj];
    et = toc;
    rt = et/j;
    eta = rt*(ns0pts-j);
    fprintf('eta: %d:%02d\r',fix(eta/60),fix(rem(eta,60)));
end
fprintf('\n');

[min_K idx] = min(K);
s0_opt = S0(idx);
plot(rl_rng,K,'-',rl_rng(idx),min_K,'r+');
title(sprintf('%s:\t tol = %2.1e',data,err_tol));
legend(sprintf('hess (k,s0): (%4.3g, %d)',rl_rng(idx(1)),min_K(1)),sprintf('proj: (%4.3g, %d)',rl_rng(idx(2)),min_K(2)));

saveas(gcf,sprintf('%s_k_vs_s0_1D',data));
% -----------------------------------------------------
    function errout = tfunc_err_2methods(Vk,Hk)
        tfunc = abs(tf_hess(Vk,Hk,r0,c,s,s0));
        err_hess = norm(tfunc_real(:) - tfunc(:),inf)/ntf;
        
        tfunc = abs(tf_proj(Vk,C,G,c,b,s));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;
        
        if err_proj < err_tol && err_hess < err_tol
            errout = [nan nan];
        else 
            errout = [err_hess err_proj] < err_tol;
        end
    end

    function errout = tfunc_err(Vk,Hk)
        tfunc = abs(tf_proj(Vk,C,G,c,b,s));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;

        if err_proj < err_tol 
            errout = nan;
        else
            errout = err_proj < err_tol;
        end
    end

end  % end main function
