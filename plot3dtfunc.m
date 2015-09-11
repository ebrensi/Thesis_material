function [h_FR h_surf h_poles h_polewt] = plot3dtfunc(data_filename,bt,ct,specs)
% [h_FR h_surf h_poles h_polewt] = plot3dtfunc(data_filename,bt,ct,specs)
%
%  plot the magnitude of the transfer function over a rectanglular region of the complex plane.
%  specs are N, projection_method, and s0


N = 136;       % ROM size
% projection_method = 'explicit';
projection_method = 'implicit';
s0 = (1+1i)*1e10;  % expansion point for the basic 1 point shift-invert Arnoldi for MOR
s0string = '(1+1i)*1e10';
nrpts = 200;   rrange = [0 10];
nipts = 200;   irange = [0 10.1];

ROI = [9 10];
FRdomain = 1i * 10.^ROI;

eigtol = 1e-12;  % any eigs (of (A-s0E)\E) with magnitude smaller than this are considered zero
reflect_about_real_axis = false;

if exist('specs','var')  
	N = specs.N;       % ROM size
	projection_method = specs.projection_method;
	s0 = specs.s0;  
end

model_name = sprintf('%ss_{%d%d}',data_filename,bt,ct); 

% define r-axis sample points
number_of_positive_rpts = floor(nrpts/2);
[~, rfrq] = tfdomain(number_of_positive_rpts,rrange);
rfrq = [fliplr(rfrq) 0 -rfrq];  % Make this plot symmetric about i-axis.
real_0 = rfrq==0;  % This is the index for r=0.


% define i-axis sample points
[~, ifrq] = tfdomain(nipts,irange,'linear');
ROImask = 10^ROI(1) <= ifrq & ifrq <= 10^ROI(2);
ROIfrq = ifrq(ROImask);


% construct the rectangular grid of real "frequency" points
[rpts ipts] = meshgrid(rfrq,ifrq);

% transfer function domain in the complex plane
s = complex(2*pi*rpts,2*pi*ipts);   % s is a nipts x nrpts matrix.


%% Construct suitably accurate ROM for analysis
[A E b c] = realization(data_filename,bt,ct);
[AE R] = make_SI_op(A,E,b,s0);
result = band_Arnoldi(N,R,AE);  % use band-Arnoldi to get projection basis V
V = result.V;
H = result.H;
rho = result.rho;

% evaluate transfer function at points s
if strcmp(projection_method,'explicit')
	if ~isreal(V)
		V = make_basis_real(V);
	end
	n = size(V,2);
	A = V'*A*V; E = V'*E*V; b = V'*b; c = V'*c;  % Replace realization with size n ROM realization
	tf = abs(transfer_function([],A,E,c,b,s ));
	tf_s0 = abs(transfer_function([],A,E,c,b,s0));
else
	n = size(V,2);
	c = V'*c;
	tf = abs(tf_implicit(H,rho,c,s0,s));
	tf_s0 = abs(tf_implicit(H,rho,c,s0,s0));
end


% plot 2d cross-section over the i-axis
h_FR = figure();
ROM_frq_response = tf(ROImask,real_0);
semilogy(ROIfrq,ROM_frq_response,'k');



% plot URM frequency response, if available
[URM_frq_response URM_ROIfrq] = URM_freq_response(data_filename,bt,ct,ROIfrq);
URM_frq_response = abs(URM_frq_response)';
err = norm(URM_frq_response - ROM_frq_response) / norm(URM_frq_response);
h = line(URM_ROIfrq,URM_frq_response,'Color','r','LineWidth',1,'LineStyle','--','Marker','+');
uistack(h,'bottom');
title(sprintf('%s Gain:  n=%g, %s,  s_0 = %s, err = %g',model_name,n,projection_method,s0string,err));
xlabel('\itf'); ylabel('|H(2\pi \itf)|');


%% plot transfer function as surface over C
log_tf = log10(tf);
h_surf = figure('name',sprintf('%s surf',model_name));
surf(rpts,ipts,log_tf);
if reflect_about_real_axis
	hold on; surf(rpts,-ipts,log_tf); hold off; % mirror
end

hold on
s0r = real(s0)/(2*pi);
s0i = imag(s0)/(2*pi);
scatter3(s0r,s0i,log10(tf_s0), 256,'r','o');           % expansion point s0
h = line([s0r s0r],[s0i s0i],zlim,'lineWidth',2,'Color','b','LineStyle','--');
uistack(h,'bottom');
hold off

% highlight frequency response
line(zeros(1,length(ifrq)),ifrq, log_tf(:,real_0),'color','w','linewidth',3);
line(zeros(1,nnz(ROImask)), ROIfrq, log10(ROM_frq_response),'color','m','linewidth',3);
if reflect_about_real_axis
	line(zeros(1,nnz(ROImask)),-ROIfrq,log10(ROM_frq_response),'color','m','linewidth',3);% mirror
end
title(sprintf('%s Gain:  n=%g, %s, s_0 = %s',model_name,n,projection_method,s0string));
ylabel('Im(s)'); xlabel('Re(s)'); zlabel('log_{10}|H(s)|');
set(gca,'XScale','linear','YScale','linear','ZScale','linear');
if reflect_about_real_axis
	centered_axis(rfrq,[-ifrq ifrq]);
else
	centered_axis(rfrq,ifrq);
end
view(68,6);
cameratoolbar;

h_poles = figure;
% Plot contour plot of full transfer function
contour(2*pi*rpts,2*pi*ipts,log10(tf));
set(gca,'XScale','linear','YScale','linear');
% if reflect_about_real_axis
% 	hold on; contour(2*pi*rpts,-2*pi*ipts,log10(tf)); hold off;
% end

%%    ************* determine and plot poles ****************
if strcmp(projection_method,'explicit')
	[AT ET Q Z V W] = qz(A,E);  % note that we have renamed V
	clear Q Z
	mu = ordeig(AT,ET);

	%% compute residues/weights
	f = sum(abs(W'*b),2);
	g = sum(abs(V'*c),2);
	residue = f(:) .* g(:);
	% This scaling is necessary for W'*En*V = I, which is an assumption for the residue
	% computation.
	vec_scaling = diag(W'*E*V);
	lambda = 1./(mu-s0);
	inf_pole = abs(lambda) < eigtol;  % Anything this far away from s0 is infinite.
	
	mu(inf_pole) = [];                % We discard infinite Ritz values for the plot
	residue(inf_pole) = [];
	lambda(inf_pole) = [];
	vec_scaling(inf_pole) = [];
	wt = abs(residue) ./ abs(vec_scaling) ./  (1+d2S(mu,FRdomain)); %abs(real(mu));  %abs(real(mu)); % Aguirre's DPI. %
else
	[W D] = eig(H); 
	lambda = diag(D);
	inf_pole = abs(lambda) <eigtol;  % Anything this far away from the origin is infinite.

	mu = s0 + 1 ./ lambda;
	
	% Compute pole mass
	f =  sum(abs(c'*W),1)';   
	g = sum(abs(W\rho),2);
	delta = abs(s0-mu) ./ (1+d2S(mu,FRdomain));
	delta(isinf(mu)) = 1;
	wt = delta .* f .* g;
	ROM_mass = sum(wt);
end

[slambda idx] = sort(lambda,'descend');   % lambda sorted by magnitude
slamu = mu(idx);   % mu sorted by lambda magnitude (convergence order)

[wt_order idx] = sort(log10(wt),'descend');
ds = dotsize(wt_order);
mu = mu(idx);

       

hold on
ninf = d2S(mu,2*pi*FRdomain) < 1e10;
scatter(real(mu(ninf)),imag(mu(ninf)),ds(ninf),ds(ninf),'*');
h = line([0 0],2*pi*10.^ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
uistack(h,'bottom');

title(sprintf('%s poles: n=%g, %s, s_0 = %s',model_name,n, projection_method, s0string));
ylabel('Im(s)'); xlabel('Re(s)');


scatter(real(s0),imag(s0),96,'r','p','MarkerFaceColor','r');           % expansion point s0

circle(real(s0), imag(s0),abs(slamu(1)-s0),'g-.');
circle(real(s0), imag(s0),abs(slamu(10)-s0),'b-.');
circle(real(s0),-imag(s0),abs(slamu(1)-s0),'g-.');
circle(real(s0),-imag(s0),abs(slamu(10)-s0),'b-.');

hold off

a = max(abs([real(mu); imag(mu)]));
xlim([-a a]);
ylim([-a a]);
putaxis

%% plot wt magnitudes with dot size
h_polewt  = figure;
scatter(1:length(wt), wt_order, ds,ds,'*');
colorbar;
title(sprintf('%s pole weight, n=%d, %s, s_0 = %s',model_name,n,projection_method,s0string));



	function ds = dotsize(orders)
		% dot size indicates order of magnitude
		sizes = [2 8 16 32 64 128];
		% 		sizes = [1 8 16 21 32 48 64 80 110 128];
		d = length(sizes);
		orders = round(orders);
		
		maxq = max(orders);
		% 		minq = min(orders);
		
		adjusted_orders = orders-maxq+d;
		adjusted_orders(adjusted_orders<1) = 1;
		ds = sizes(adjusted_orders);
	end

	function centered_axis(rpts,ipts)
		rmax = max(rpts); rmin = min(rpts);
		imax = max(ipts); imin = min(ipts);
		rlength = rmax-rmin; ilength = imax-imin;
		
		if rlength > ilength
			rlim = [rmin rmax];
			sum = imin + imax;
			ilim = [sum-rlength  sum+rlength]/2;
		else
			ilim = [imin imax];
			sum = rmin+rmax;
			rlim = [sum-ilength  sum+ilength]/2;
		end
		xlim(rlim);
		ylim(ilim);
	end
end
