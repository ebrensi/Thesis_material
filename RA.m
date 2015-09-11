function [V H] = RA(A,E,b,S,J)
% [V H] = RA(A,E,b,S,J)
%
%  Given A,E in R^{N x N}
%        S = [s1 s2 ... sn]
%        J = [j1 j2 ... jn]
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
q = sum(J);     % total model size

V = zeros(N,q);
H = zeros(q,q);

%% Initialize
[Hi rk] = make_SI_op(A,E,b,S(1));

%% Begin RA Iterations
for i = 1:n
	for j = 1:J(i)
		k = (i-1)*J(i) + j;
		%  Generate New Orthormal Vector vk
		norm_rk = norm(rk);
		V(:,k) = rk / norm_rk;
		if k > 1
			H(k,k-1) = norm_rk;
		end
		
		
		% Update the Residue rk for the next iteration
		if j == J(i) &&  i < n
			[Hi rk] = make_SI_op(A,E,b,S(i+1));
		else
			rk = Hi(V(:,k));
		end
		
		% make rk orthogonal to all previous v's (1:k) via MGS
		for t = 1:k
			vt = V(:,t);
			H(t,k) = vt' * rk;
			rk = rk - H(t,k)*vt;
		end
		%% optional reorthogonalization 
		if reorthogonalize
			h2 = zeros(k,1);
			for t = 1:k
				vt = V(:,t);
				h2(t) = vt' * rk;
				rk = rk - h2(t)*vt;
			end
			H(1:k,k) = H(1:k,k) + h2;
		end
	end
	
end
