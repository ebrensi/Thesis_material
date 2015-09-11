function [data rng n_cvrgd swt]  = poles_vs_s0_1D(data,kmax,npts,ritz_tol)
% 

rng = linspace(1e-6,2,npts);
S0 = rng *sqrt(-1)*pi * 1e10;

[C G c b] = inputdata(data);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));
ntf = norm(tfunc_real,inf);
swt = zeros(1,npts);
n_cvrgd = zeros(1,npts);
data = zeros(kmax*npts,4);

k=0;
for j = 1:npts
    s0 = S0(j);
    [A r0] = makeA(C,G,b,s0);
    V  = arnoldiR(A,r0,kmax);
    [mu Z rr tf_terms] = approx_poles_proj(V,C,G,c,b,s0);
    cvrgd = rr < ritz_tol;
    n_cvrgd(j) = nnz(cvrgd);
    wt = tf_terms_weight(tf_terms(cvrgd,:));
    swt(j) = sum(wt);
    idx = k+(1:n_cvrgd(j));
    data(idx,1:2) = [mu(cvrgd) wt];
    data(idx,3) = j;
    data(idx,4) = rng(j);
    
    if ~isempty(idx)
        k = idx(end);
    end

    fprintf('%d ',j);
end
fprintf('\n');
data(k+1:end,:) = [];
data2 = sortrows(data,1);

% fprintf('\n');
% [max_npoles idx_maxnp] = max(n_cvrgd)
% s0_opt_proj = S0(idx_maxnp)

plot(rng,swt);
title('wt');
figure;
plot(rng,n_cvrgd);



% semilogy(rl_rng,err_proj,rl_rng(idx_proj),err_proj(idx_proj),'r+');
% title(sprintf('%s:\t m = %d    eh = %4.3g,    ep = %4.3g',data,kmax,min_err_proj));
% legend(sprintf('proj best s0: %4.3g',rl_rng(idx_proj)));
% 
% saveas(gcf,sprintf('%s_1D',data));
% % exportfig(gcf,sprintf('%s_k%d',data,kmax));

%% -------------------------------------------------------------

end
