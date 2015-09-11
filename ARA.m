function V  = ARA(A,E,b,S,J)

N = length(b);  % order of realization
n = length(S);  % number of runs "restarts"
q = sum(J);     % total model size

V = zeros(N,q);
H = zeros(q,q);
Q = [];

%% Initialize
for i = 1:n
	[Hi rk] = make_SI_op(A,E,b,S(i));
	[Vi vi_next projected_Hi op_norm_est] = augmented_arnoldi(Q,Hi,rk,J(i));
	
	[W D] = eig(projected_Hi);
	Ritz = diag(D);
	
	rr = norm(vi_next) * abs(W(end,:)).' ./ abs(Ritz); % relative residual 2-norm
	chk = sortrows([Ritz abs(Ritz) rr],2);
	
	converged = rr <= 1e-2;
	Y = Vi*W(:,converged);
	Q = qr([real(Y) imag(Y) Q]);
end

	
end
