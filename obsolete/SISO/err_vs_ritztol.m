function [rt err rnk] = err_vs_ritztol(data,m,K,npts)

rt = logspace(-16,1,npts);
good = [0.4300  0.6600  0.8629  1.2500 ];

err = zeros(npts,1);
rnk = zeros(npts,1);
for j = 1:npts
    j
    [err(j) V] = AT4(data,m,K,rt(j),good);
    rnk(j) = rank(V);
end
loglog(rt,err);
figure;
semilogx(rt,rnk);
