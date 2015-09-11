function cmpr(datafilename,nvals,s0vals,nsamps)

markers = {'+','o','.','s','d','p','h','*','x'};
s0trials = [pi*1e10 (1+1i)*pi*1e10 1000+1i*pi*1e10];
names = {'\pi10^{10}','\pii(1+1i)10^{10}','\pii10^{10}+10^3','multi','multi2',...
	'multi3','multi4','multi5'};

nmax = sum(nvals(1,:));
nfirst = 2;

if ~exist('nsamps','var') || isempty(nsamps) || nsamps > nmax-nfirst +1
	nsamps = nmax-nfirst+1;
end

numberofmodels = size(nvals,1);
cplx_errs_rstrt = zeros(nsamps,numberofmodels);
cplx_ROMsize_rstrt = zeros(nsamps,numberofmodels);

for j=1:numberofmodels
	result = nrmods(datafilename,nvals(j,:),s0vals(j,:),nsamps);
	s0vals(j,:) = result.s0vals;
	cplx_errs_rstrt(:,j) = result.errs(:,1);
	cplx_ROMsize_rstrt(:,j) = result.cplx_ROMsize;
	nmax = max(cplx_ROMsize_rstrt);
end


%% now we compute several single-run ROMs with various s0 values.
n_s0 = numel(s0trials);
cplx_errs_fixed = zeros(nsamps,n_s0);
cplx_ROMsize_fixed = zeros(nsamps,n_s0);

for j = 1:n_s0
	s0 = s0trials(j);
	result = nrmods(datafilename,nmax,s0,nsamps);
	cplx_errs_fixed(:,j) = result.errs(:,1);
	cplx_ROMsize_fixed(:,j) = result.cplx_ROMsize;
end

figure('name','cplx prj errs');
h = semilogy(cplx_ROMsize_fixed,cplx_errs_fixed,cplx_ROMsize_rstrt,cplx_errs_rstrt);
set(h,{'Marker'},markers(1:n_s0+numberofmodels)');
legend(h,names(1:n_s0+numberofmodels));
title('tfunc errors (cplx proj) vs ROM dimension')

