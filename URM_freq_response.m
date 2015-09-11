function [FR frq s] = URM_freq_response(data_filename,bt,ct,frq)
% [FR frq s] = URM_freq_response(data_filename,bt,ct,frq)
%
% compute the frequency response of the undreduced model for input values of frq
% via the original formulation
%  H(s) = C' * (sE-A)^-1 * B

if ~exist('frq','var') || isempty(frq)
	[s frq] = tfdomain();
end

 SISO = exist('bt','var') & ~isempty(bt) & exist('ct','var') & ~isempty(ct);

npts = length(frq);
f1 = round(log10(frq(1)));
f2 = round(log10(frq(end)));
s = 2i*pi*frq;

if SISO
	model_name = sprintf('%ss%d%d',data_filename,bt,ct);
else
	model_name = data_filename;
end

FR_filename = sprintf('FR_%s_%d_%d_%d.mat',model_name,npts,f1,f2);


if exist(FR_filename,'file') == 2
	vars = load(FR_filename);
	FR = vars.FR;
else
	if SISO
		[A E C B] = realization(data_filename,bt,ct);
	else
		[A E C B] = realization(data_filename);
	end
	FR = transfer_function([],A,E,C,B,s);
	save(FR_filename,'FR','frq');
end
