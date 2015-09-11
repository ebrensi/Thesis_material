function [mu wt ROM] = modal_trunc(data,num_modes)
% [mu wt ROM] = modal_trunc(data,num_modes)
%
% Pole/eigenvalue analysis of the unreduced model transfer function
%  H(s) = c'(sE-A)\b

animation_on = false;

%% ** Read in realization for this model
[A E c b] = realization(data);
[s frq] = tfdomain();


%% determine poles of this model, which are eigenvalues of (A,E)
[AA,EE,Q,R,Z,W] = examp_qz(data);
mu = ordeig(AA,EE);

% Just in case there are some indeterminant eigs
ill_disposed = isnan(mu);
if any(ill_disposed)
	mu(ill_disposed) = [];
	Z(:,ill_disposed) = [];
	W(:,ill_disposed) = [];
	warning('there were ill disposed eigs');
end

% We scale eigenvectors W and V so that W'*E*V = I.
eig_vec_scaling = abs(diag(W'*E*Z));

% Compute residues
f = Z'*c;
g = W'*b;
residue = conj(f) .* g ./ eig_vec_scaling;


%  Compute pole weights
infinite_poles = isinf(mu);
% wt(infinite_poles) = abs(residue(infinite_poles));
wt(infinite_poles) = 0;
wt(~infinite_poles) = abs(residue(~infinite_poles)) ./ d2S(mu(~infinite_poles),s);

% sort poles by descending weight
[wt_sorted dom_pole_idx] = sort(wt,'descend');
 
% Now we construct a ROM by projecting on to the invariant subspace associated with the
%  most dominant poles.
if ~exist('num_modes','var')
	num_modes = length(mu);
end
dpoles = dom_pole_idx(1:num_modes);

figure;
 plot_poles(mu(dpoles), wt(dpoles),s)


%  Re-order the schur decomposition to have dominant poles in the upper submatrix. 
mask = false(length(mu),1);
mask(dom_pole_idx(1:num_modes)) = true;
num_modes = nnz(mask);
[AAS,EES,QS,ZS] = ordqz(AA,EE,Q,R,mask);

% Oblique projection of URM onto shur space of dominant poles. 
n = num_modes;
cp = c'*ZS(:,1:n);
bp = QS(1:n,:)*b;
Ap = AAS(1:n,1:n);
Ep = EES(1:n,1:n);

%  Construct projected model transfer function from projected realization.
Ls = length(s);
ROM_tfunc = zeros(1,Ls);
optU.UT = true; 
for j = 1:Ls
    ROM_tfunc(j) = cp*linsolve(s(j)*Ep-Ap,bp,optU);
end


% compute the urm transfer function explicitly via standard formulation
URM_tfunc = examp_tfunc(data,frq);

ntf = norm(URM_tfunc,inf);
err = norm(ROM_tfunc - URM_tfunc,inf)/ntf;

% plot both of them
figure;
loglog(frq,abs(URM_tfunc),'r',frq,abs(ROM_tfunc));
set(gca,'Xscale', 'linear');
title(sprintf('ModalTrunc %s Gain:  n=%g,     err: %g,      mesh:%d',data,n,err,length(frq)));

% figure;
% loglog(frq,abs(URM_tfunc - ROM_tfunc)/norm(URM_tfunc));
% title(sprintf('Diff   num modes: %g,  rel err: %g',n,err));


if animation_on
	figure
	
	% *** plot transfer function using successively more modes ****
	err = zeros(length(wt),1);
	done = false;
	j = 1;
	while ~done && j < length(wt)
		Hj = abs(sum(tf_terms(1:j,:),1));
		err(j) = norm(Hj - tfunc_exp,inf)/ntf;
		loglog(frq,abs(tfunc_exp),'r',frq,Hj);
		title(sprintf('1:%d, err = %g',j,err(j)));
		drawnow;
		done = waitforbuttonpress;
		j = j+1;
	end
end


%% -----------------------------------------------------------------------
	function [mu wt] = poles_via_shift_invert(C,G,c,b,s0)
		%% ****** compute eigenvalues of H(s0) = (C-s0*G)^-1 * G  ******
		s0G_C = s0*G-C;
		r0 = s0G_C\b;
		% ** compute eigs of H **
		H = -s0G_C\G;
		[Z L] = eig(full(H));
		L = diag(L);
		mu = s0 + 1./L;
		
		%% compute transfer function terms
		f = Z'*c;
		g = Z\r0;
		
		fg = conj(f) .* g;
		% 		tf_terms = bsxfun(@rdivide, fg, 1-L*(s-s0));
		
		wt = abs(fg) .* abs(s0-mu) ./ d2S(mu);
	end
end % main function
