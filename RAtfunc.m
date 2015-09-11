function RAtfunc(data_filename,S,J)
% tf = RAtfunc(data_filename,S,J)
%
% computes the ROM transfer function of a model using Rational Arnoldi (RA) method 

% Get the domain for the transfer function 
[domain frq] = tfdomain();
nplotpts = length(domain);

% Input data (matrices C,G,b) from data file; note H(s) = c'*(C-s*G)*b
[A E c b] = realization(data_filename);

% Construct rational Krylov projection basis using RA
[V H] = RA(A,E,b,S,J);


% realify basis
if ~isreal(V)
	V = make_basis_real(V);
end
n = size(V,2);

% Explicitly compute orthogonally projected ROM
% Evaluate the transfer function over the given domain
ROM_tfunc = transfer_function(V,A,E,c,b,domain);

% Evaluate URM transfer function over the given domain
URM_tfunc = examp_tfunc(data_filename,frq);

err = norm(ROM_tfunc - URM_tfunc) / norm (URM_tfunc); 
% err_max = max(abs(ROM_tfunc - URM_tfunc)) / max(abs(URM_tfunc)); 

% plot both of them
figure;
loglog(frq,abs(URM_tfunc),'r',frq,abs(ROM_tfunc));
title(sprintf('RA %s Gain:  n=%g,     err: %g,      mesh:%d',data_filename,n,err,nplotpts));

% figure;
% loglog(frq,abs(URM_tfunc - ROM_tfunc)/norm(URM_tfunc));
% title(sprintf('%s diff:  n=%g,     err: %g,      mesh:%d',data_filename,n,err,nplotpts));

