function s0 = AORAtfunc(data_filename,S,q)
% tf = AORAtfunc(data_filename,S,q)
%
% computes the ROM transfer function of a model using Rational Arnoldi (RA) method 


% Get the domain for the transfer function 
[domain frq] = tfdomain();
nplotpts = length(domain);

% Input data (matrices C,G,b) from data file; note H(s) = c'*(C-s*G)*b
[A E b c] = realization(data_filename);

% Construct rational Krylov projection basis using RA
[V s0] = AORA(A,E,c,b,S,q);

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
% set(gca,'Xscale', 'linear');
% ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
% 		   'YAxisLocation','right',...
%            'Color','none',...
% 		   'XColor', 'red',...
% 		  'XLim', log10((10.^freq_interval)*1e-9),...
%            'Xscale','linear',...
% 		   'XMinorTick','on',...
% 	   'XLimMode','manual');
title(sprintf('AORA  %s Gain:  n=%g,     err: %g,      mesh:%d',data_filename,n,err,nplotpts));

% figure;
% loglog(frq,abs(URM_tfunc - ROM_tfunc)/norm(URM_tfunc));
% title(sprintf('%s diff:  n=%g,     err: %g,      mesh:%d',data_filename,n,err,nplotpts));
