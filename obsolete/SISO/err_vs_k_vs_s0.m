function [err_h err_p] = err_vs_k_vs_s0(data,nkpts,ns0pts)
% [err_h err_p] = err_vs_k_vs_s0(data,nkpts,ns0pts)

kmax = 200;
[C G c b] = inputdata(data);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);

rl_rng = linspace(1e-6,2,ns0pts);
S0 = rl_rng *pi * 1e10;

kvals = fix(linspace(2,kmax,nkpts));
err_h = zeros(nkpts,ns0pts);
err_p = zeros(nkpts,ns0pts);

tic;
for j = 1:ns0pts
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    err = log10(step_through_arnoldi(A,r0,kvals,@tfunc_err,2));
    err_h(:,j) = err(:,1);
    err_p(:,j) = err(:,2);
    
    et = toc;
    rt = et/j;
    eta = rt*(ns0pts-j);
    fprintf('eta: %d:%02d\r',fix(eta/60),fix(rem(eta,60)));
end

[X Y] = meshgrid(rl_rng,kvals);

% plot err and opt_err for Hessenberg formulation
[min_y iy] = min(err_h);
[min_err ix] = min(min_y);
xmin = ix; ymin = iy(ix);
err_opt_h = err_h(ymin,xmin);
opts0_h = rl_rng(xmin);
optk_h = kvals(ymin);

mesh(X,Y,err_h);
hold on
plot3(opts0_h,optk_h,err_opt_h,'r+');
hold off
title(sprintf('%s (hess): log(err) = %4.3g at \t (s0,k) = (%2.3g, %d)',data,err_opt_h,opts0_h,optk_h));
xlabel('s0'); ylabel('k'); zlabel('log err');

figure;

% plot err and opt_err for PRIMA formulation
[min_y iy] = min(err_p);
[min_err ix] = min(min_y);
xmin = ix; ymin = iy(ix);
err_opt_p = err_p(ymin,xmin);
opts0_p = rl_rng(xmin);
optk_p = kvals(ymin);

mesh(X,Y,err_p);
hold on
plot3(opts0_p,optk_p,err_opt_p,'r+');
hold off
title(sprintf('%s (PRIMA): log(err) = %4.3g at \t (s0,k) = (%2.3g, %d)',data,err_opt_p,opts0_p,optk_p));
xlabel('s0'); ylabel('k'); zlabel('log err');

%% ---------------------------------------------------------
    function err = tfunc_err(Vk,Hk)
        tfunc = abs(tf_hess(Vk,Hk,r0,c,s0));
        err_hess = norm(tfunc_real(:) - tfunc(:),inf)/ntf;
        
        tfunc = abs(tf_proj(Vk,C,G,c,b));
        err_proj = norm(tfunc_real(:) - tfunc(:))/ntf;
        
        err = [err_hess err_proj];
    end
end
