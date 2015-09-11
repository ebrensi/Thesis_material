function k_vs_s0_2D(data,nkpts,ns0pts)
% determines number of steps k<kmax required to achieve a tfunc error of err_tol
%  for several values of s0

err_tol = 0.01;
kmax = 200;
kvals = fix(linspace(2,kmax,nkpts));

rl_rng = linspace(1e-6,2,ns0pts);
im_rng = linspace(0,2,ns0pts);

[X Y] = meshgrid(rl_rng,im_rng);
S0 = (X + Y*sqrt(-1)) *pi * 1e10;

[C G c b] = inputdata(data);
frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;

tfunc_real = abs(tfunc_urm(data,C,G,c,b))';
ntf = norm(tfunc_real,inf);

K_hess = zeros(size(S0));
K_proj = zeros(size(S0));
nel = length(S0(:));
tic;

for j = 1:nel
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    err = step_through_arnoldi(A,r0,kvals,@tfunc_err,2);
    k_hess = find(err(:,1),1);
    if isempty(k_hess) 
        k_hess = kmax;
    end
    k_proj = find(err(:,2),1);
    if isempty(k_proj) 
        k_proj = kmax;
    end
    K_hess(j) = k_hess; 
    K_proj(j) = k_proj;
    et = toc;
    rt = et/j;
    eta = rt*(nel-j);
    fprintf('eta: %d:%02d\r',fix(eta/60),fix(rem(eta,60)));
end
fprintf('\n');

[min_y iy] = min(K_hess);
[min_K_hess ix] = min(min_y);
xmin = ix; ymin = iy(ix);
s0_opt_hess = S0(ymin,xmin);
optx = rl_rng(xmin);
opty = im_rng(ymin);
mesh(X,Y,K_hess);
hold on
plot3(optx,opty,min_K_hess,'r+');
hold off
title(sprintf('%s(hess): opt s0 = (%3.2g, %3.2gi)  k = %d, \t err tol = %2.1e',data,optx,opty,min_K_hess,err_tol));
xlabel('Re'); ylabel('Im');
saveas(gcf,sprintf('%s_%s_2D',data,'hess'));

[min_y iy] = min(K_proj);
[min_K_proj ix] = min(min_y);
xmin = ix; ymin = iy(ix);
s0_opt_proj = S0(ymin,xmin);
optx = rl_rng(xmin);
opty = im_rng(ymin);
mesh(X,Y,K_proj);
hold on
plot3(optx,opty,min_K_proj,'r+');
hold off
title(sprintf('%s(PRIMA): opt s0 = (%3.2g, %3.2gi)  k = %d, \t err tol = %2.1e',data,optx,opty,min_K_hess,err_tol));
xlabel('Re'); ylabel('Im');
saveas(gcf,sprintf('%s_%s_2D',data,'proj'));

% -----------------------------------------------------
    function errout = tfunc_err(Vk,Hk)
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

end  % end main function
