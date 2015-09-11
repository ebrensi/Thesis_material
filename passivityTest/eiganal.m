function [ee tt] = eiganal(data)
load(sprintf('poles_%s',data));
ee = [sort(mu_CG) sort(1./mu_GC) sort(mu_ST) sort(1./mu_TS)];
lbl = {'(C,G)' '(G,C)' '(S,T)' '(T,S)'};
tt = zeros(size(ee));

for i=1:4
    mu = process_eigs(ee(:,i));
    tt(1:length(mu),i) = mu;
%     figure;
%     hist(mu,20);
%     title(sprintf('%s: %s',data,char(lbl(i))));
end

last_nonzero = 1+find(~any(tt(2:end,:),2),1,'first');
tt(last_nonzero:end,:) = [];

semilogy(abs(tt(1:last_nonzero-1,:)));
legend(lbl);
