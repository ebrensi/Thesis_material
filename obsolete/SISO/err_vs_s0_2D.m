function s0_opt = err_vs_s0_2D(data,kmax,npts,meth)
% generates reduced models for several expansion points s0 to see which
% work best.

rl_rng = linspace(1e-6,2,npts);
im_rng = linspace(0,2,npts);

[X Y] = meshgrid(rl_rng,im_rng);
S0 = (X + Y*sqrt(-1)) *pi * 1e10;

[C G c b] = inputdata(data);

tfunc_real = abs(tfunc_urm(data,C,G,c,b))';
ntf = norm(tfunc_real,inf);

err = zeros(size(S0));
nel = length(S0(:));
tic;
for j = 1:nel
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    [V H] = arnoldi(A,r0,kmax);
    if strcmp(meth,'hess')
        tfunc = abs(tf_hess(V,H,r0,c,s0));
    else
        tfunc = tf_proj(V,C,G,c,b);
    end
    
    err(j) = norm(tfunc(:)-tfunc_real(:))/ntf;
    if mod(j,10) == 0
        et = toc;
        rt = et/j;
        eta = rt*(nel-j);
        fprintf('eta: %d:%d\r',fix(eta/60),fix(rem(eta,60)));
    end
end
fprintf('\n');

[min_y iy] = min(err);
[min_err ix] = min(min_y);
xmin = ix; ymin = iy(ix);
s0_opt = S0(ymin,xmin);
err_opt = err(ymin,xmin);
optx = rl_rng(xmin);
opty = im_rng(ymin);

mesh(X,Y,log10(err));
hold on
plot3(optx,opty,log10(err_opt),'r+');
hold off

title(sprintf('%s(%s): opt s0 = (%3.2g, %3.2gi)    err = %3.2e    k = %d',data,meth,optx,opty,err_opt,kmax))
xlabel('Re'); ylabel('Im');
saveas(gcf,sprintf('%s_%s_2D',data,meth));

%% -------------------------------------------------------------

end
