function [err_proj] = err_vs_s0_1D(data,kmax,npts)
% generates reduced models for several expansion points s0 to see which
% work best.

rl_rng = linspace(1e-6,2.5,npts);
S0 = rl_rng *sqrt(-1)*pi * 1e10;

[C G c b] = inputdata(data);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);
err_proj = zeros(1,npts);

for j = 1:npts
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    [V H] = arnoldi(A,r0,kmax);
    tfunc_proj = tf_proj(V,C,G,c,b);

    fprintf('%d ',j);
    err_proj(j) = norm(tfunc_proj-tfunc_real,inf)/ntf;
end
fprintf('\n');
[min_err_proj idx_proj] = min(err_proj);
s0_opt_proj = S0(idx_proj);

semilogy(rl_rng,err_proj,rl_rng(idx_proj),err_proj(idx_proj),'r+');
title(sprintf('%s:\t m = %d    eh = %4.3g,    ep = %4.3g',data,kmax,min_err_proj));
legend(sprintf('proj best s0: %4.3g',rl_rng(idx_proj)));

saveas(gcf,sprintf('%s_1D',data));
% exportfig(gcf,sprintf('%s_k%d',data,kmax));

%% -------------------------------------------------------------

end
