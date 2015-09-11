function [V s0] = AORA(A,E,c,b,S,q)
% [V s0] = AORA(A,E,c,b,S,q)
%
%  Given A,E in R^{N x N}
%        S = [s1 s2 ... sn]
%
%  This is the Rational Arnoldi process for computing a basis for the union of n Krylov
%  subspaces K_ji(Hi,ri) for i = 1...n
%
%      where Hi = (A-si E)\E  and ri = (si E - A)\b
%
%  i.e. we perform ji iterations with each expansion point si
%
%
% Output:  n x n matrix H and basis vectors v as columns of N x q matrix V.
%           where q = sum(J)
%
% This was adapted from pseudocode from Table 1 of
%  "An adaptive Arnoldi method for model-order reductions of linear time-invariant systems"
%   (2004) by Lee, Chu, and Feng
%
% Efrem Rensi 8/2013

reorthogonalize = false;

N = length(b);  % order of realization
n = length(S);  % number of runs "restarts"

V = zeros(N,q);
H = zeros(q,q,n);
r = zeros(N,n);
k = zeros(N,n);
s0 = zeros(q,1);
hpi = zeros(1,n);

%% Initialize
for i = 1:n
	[op{i} r(:,i)] = make_SI_op(A,E,b,S(i));
	k(:,i) = r(:,i);
	hpi(i) = 1;
end

%% Begin AORA Iterations
for j = 1:q
	% Select the Expansion point with the Maximum Output Moment Error
	moment_err = abs(c'*r .* hpi);
	[~, idx_maxerr] = max(moment_err);
	s0(j) = S(idx_maxerr);
	
	% Generate Orthonormal Vector at s0(j)
	nrm_r = norm(r(:,idx_maxerr));
	if j > 1
		H(j,j-1,idx_maxerr) = nrm_r;
	end
	V(:,j) = r(:,idx_maxerr) / nrm_r;
	hpi(idx_maxerr) = hpi(idx_maxerr) * nrm_r;
	
	% Update the Residue for the next iteration
	for i = 1:n
		if i == idx_maxerr
			k(:,i) = op{i}(V(:,j));
		end
		r(:,i) = k(:,i);
		for t = 1:j
			vt = V(:,t);
			H(t,j,i) = vt' * r(:,i);
			r(:,i) = r(:,i) - H(t,j,i)*vt;
		end
		%% optional reorthogonalization 
		if reorthogonalize
			h2 = zeros(j,1);
			for t = 1:j
				vt = V(:,t);
				h2(t) = vt' * r(:,i);
				r(:,i) = r(:,i) - h2(t)*vt;
			end
			H(1:j,j,i) = H(1:j,j,i) + h2;
		end
	end
end
