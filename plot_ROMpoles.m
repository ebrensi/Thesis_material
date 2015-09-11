function h = plot_ROMpoles(mu,s0,rr,mass,ROI,keep)
cax = [-9 0];
ds_mass = dotsize(mass);
dc_rr = log10(rr);  
dc_rr(dc_rr < -16) = -16;
dc_rr(dc_rr > 0) = 0;

if any(isnan(dc_rr))
	dbstop;
end

%% 	********** plot poles of the transfer function sized by wt ******************
scatter(real(mu),imag(mu),ds_mass,dc_rr,'*');
caxis(cax);

hold on
if nnz(keep) > 0
	scatter(real(mu(keep)),imag(mu(keep)),ds_mass(keep),'k','o');
end
scatter(real(s0),imag(s0),128,'k','x');    % interpolation point s0
for j = 1:length(s0)
	[dist_to_closest idx_of_closest] = min(abs(s0(j) - mu));
	mu_closest = mu(idx_of_closest);
	circle(real(s0(j)),imag(s0(j)), dist_to_closest ,'b-.');
end
hold off
h = line([0 0],2*pi*10.^ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
uistack(h,'bottom');
% axis(gca,ax);
% putaxis


	function ds = dotsize(mass)
		orders = log10(mass);
		% dot size indicates order of magnitude
		sizes = [8 16 32 64 128];
		d = length(sizes);
		orders = round(orders);
		maxq = max(orders);		
		adjusted_orders = orders-maxq+d;
		adjusted_orders(adjusted_orders<1) = 1;
		ds = sizes(adjusted_orders);
	end

	function dc = dotcolor(rr)
		% maps log(rr) in [-1e9 1]   to   dc in [0 1] 
		dc = log10(rr);
		dc(dc<-9) = -9;
		dc(dc>1) = 1;
	end
end
