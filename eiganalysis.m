function eiganalysis(datafile,bt,ct)

s0 = (.1 + 5i)*2*pi*1e9;  % expansion point for the basic 1 point shift-invert Arnoldi for MOR
s0string = '(.1 + 5i)*2*\pi*1e9';

N = 100;       % ROM size
projection_method = 'explicit';
% projection_method = 'implicit';
ROI = [9 10];
FRdomain = 10.^ROI*1i;

eigtol = 1e-12;  % any eigs with magnitude smaller than this are considered zero.

model_name = sprintf('%ss_{%d%d}',datafile,bt,ct); 
[A E C B] = realization(datafile,bt,ct);
[AE R] = make_SI_op(A,E,B,s0);

result = band_Arnoldi(N,R,AE);
V = result.V;
H = result.H;
rho = result.rho;


if strcmp(projection_method,'explicit')
	if ~isreal(V)
		V = make_basis_real(V);
	end
	n = size(V,2);
	A = V'*A*V; E = V'*E*V; b = V'*B; c = V'*C;  % Replace realization with size n ROM realization

	[AT ET Q Z V W] = qz(A,E);  % note that we have renamed V
	clear Q Z
	mu = ordeig(AT,ET);
	
	%% compute residues/weights
	f = sum(abs(W'*b),1);
	g = sum(abs(V'*c),1);
	residue = conj(f) .* g;
	
	% This scaling is necessary for W'*En*V = I, which is an assumption for the residue
	% computation.
	vec_scaling = diag(W'*E*V);
	lambda = 1./(mu-s0);
	inf_pole = abs(lambda) < eigtol;  % Anything this far away from s0 is infinite.
	mu(inf_pole) = [];                % We discard infinite Ritz values for the plot
	residue(inf_pole) = [];
	lambda(inf_pole) = [];
	vec_scaling(inf_pole) = [];
	mass = abs(residue) ./ abs(vec_scaling) ./  (d2S(mu,FRdomain)+1); %abs(real(mu)); % Aguirre's DPI. %
else
	[W D] = eig(H); 
	n = length(W);
	
	lambda = diag(D);
	inf_pole = abs(lambda) < eigtol;  % Anything this far away from the origin is infinite.
	lambda(inf_pole) = [];
	W(:,inf_pole) = [];
	
	mu = s0 + 1 ./ lambda;
	f =  sum(abs(C'*V*W),1)';   
	g = sum(abs(W\rho),2);
	delta = abs(s0-mu) ./ d2S(mu,FRdomain);
	mass = delta .* f .* g;
end

[slambda idx] = sort(lambda,'descend');   % lambda sorted by magnitude
slamu = mu(idx);                          % mu sorted by magnitude of lambda

order_of_lambda = log10(abs(slambda));

%% plot eigenvalues lambda of (A-s0E)\E
h_lambda = figure;
scatter(real(slambda),imag(slambda),32,order_of_lambda,'*');
tit = 'Ritz \lambda of (A-s_0E)^{-1}E';
title(sprintf('%s: %s, n=%d, %s,  s_0 = %s',model_name,tit,n,projection_method,s0string));


%% plot lambda order of magnitude with dot size
h_eigmag  = figure;
scatter(1:length(slambda), order_of_lambda, 32,order_of_lambda,'*');
colorbar;
tit = 'log_{10}|\lambda|';
title(sprintf('%s: %s, n=%d, %s, s_0 = %s',model_name,tit,n,projection_method,s0string));




% plot eigenvalues mu of (A,E), with s0 and indicate which will be the first to converge
%  i.e.  those associated with large lambda
h_mu = figure;
hold on
scatter(real(slamu),imag(slamu),32,order_of_lambda,'*');  % all finite Ritz-values
scatter(real(s0),imag(s0),256,'r','x');   % expansion point (shift) s0
circle(real(s0),imag(s0),abs(slamu(1)-s0),'g-.');
circle(real(s0),imag(s0),abs(slamu(10)-s0),'g-.');
hold off
h = line([0 0],2*pi*10.^ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
uistack(h,'bottom');
a = max(abs([real(mu); imag(mu)]));
xlim([-a a]); ylim([-a a]);
putaxis
tit = 'Ritz \mu of (A,E)';
title(sprintf('%s: %s, n=%d with s_0 = %s',model_name,tit,n,s0string));
ax = axis(gca);

if size(B,2) == 1 && size(C,2) == 1;  % is this a SISO model?	
	%% plot wt magnitudes with dot size
h_polewt  = figure;
wt_order = log10(abs(mass));
[swt_order idx] = sort(wt_order,'descend');

ds = dotsize(swt_order);
scatter(1:length(mass), swt_order, ds,ds,'*');
colorbar;
title(sprintf('%s pole weight, n=%d, %s',model_name,n,projection_method));
	
	
%% 	********** plot poles of the transfer function sized by wt ******************
figure;
hold on;
line([0 0],2*pi*ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
scatter(real(mu(idx)),imag(mu(idx)),ds,ds,'*'); % poles by size according to relative order
scatter(real(s0),imag(s0),96,'r','p','MarkerFaceColor','r');           % expansion point s0
circle(real(s0),imag(s0), abs(slamu(1)-s0),'g-.');
circle(real(s0),imag(s0), abs(slamu(10)-s0),'b-.');
circle(real(s0),-imag(s0),abs(slamu(1)-s0),'g-.');
circle(real(s0),-imag(s0),abs(slamu(10)-s0),'b-.');

plot(0,2*pi*ROI(1),'k+', 0, 2*pi*ROI(2),'k+');   % delimiting makers for segment
tit = 'Ritz-poles \mu of C^T(A-sE)^{-1}B';
title(sprintf('%s:  %s, %s,  n=%d,   s_0 = %s',model_name,tit,projection_method,n,s0string));
hold off
h = line([0 0],2*pi*10.^ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
uistack(h,'bottom');
axis(gca,ax);
putaxis	
end

	function ds = dotsize(orders)
		% dot size indicates order of magnitude
		sizes = [2 8 16 32 64 128];
%         sizes = [1 8 16 21 32 48 64 80 110 128];
		d = length(sizes);
		orders = round(orders);
		
		
		maxq = max(orders);
		% 		minq = min(orders);
		
		adjusted_orders = orders-maxq+d;
		adjusted_orders(adjusted_orders<1) = 1;
		ds = sizes(adjusted_orders);
	end

end % main function

