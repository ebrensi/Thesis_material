function abort = tfunc_surf_exp(rzn,s0,gridsize,ax,ROI)
%  tfunc_surf_exp(A,E,B,C,s0)
%
%  Plot the magnitude of the transfer function over a rectanglular rgion of the complex plane.
%  specs are N, projection_method, and s0

A = rzn.A;
E = rzn.E;
B = rzn.B;
C = rzn.C;

N = size(A,2);
nrpts = gridsize(1);   
nipts = gridsize(2);   
rrange = ax(1:2); 
irange = ax(3:4);

FRdomain = 1i * 10.^ROI;



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


h = waitbar(0,'constructing surface plot...');
abort = false;
% evaluate transfer function at points s
for i = 1:size(s,1)
	tf(i,:) = abs(transfer_function([],A,E,C,B,s(i,:)));
	if ~ishandle(h)
		abort = true;
		break;
	else
		waitbar(i/size(s,1),h);
	end
end

if abort
	return;
end
delete(h);
tf_s0 = abs(transfer_function([],A,E,C,B,s0));
ROM_frq_response = tf(ROImask,real_0);


%% plot transfer function as surface over C
log_tf = log10(tf);
h_surf = figure;
surf(rpts,ipts,log_tf);

hold on
s0r = real(s0)/(2*pi);
s0i = imag(s0)/(2*pi);
h = scatter3(s0r,s0i,log10(tf_s0), 1024,'g','*');           % expansion point s0
uistack(h,'top');

for i = 1:length(s0)
	h = line([s0r(i) s0r(i)],[s0i(i) s0i(i)],zlim,'lineWidth',2,'Color','b','LineStyle','--');
end
hold off

% highlight frequency response
line(zeros(1,length(ifrq)),ifrq, log_tf(:,real_0),'color','w','linewidth',3);
line(zeros(1,nnz(ROImask)), ROIfrq, log10(ROM_frq_response),'color','m','linewidth',3);
set(gca,'XScale','linear','YScale','linear','ZScale','linear');
centered_axis(rfrq,ifrq);
view(68,6);
cameratoolbar;


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
