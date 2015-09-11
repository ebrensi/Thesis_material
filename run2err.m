function ROM_tfunc = run2err(data_filename,s0,tfunc_tol)

% parameters
dtol = sqrt(eps);
ctol = dtol;
ROI = [8 10];
stepsize = 10;

% Model setup
model_name = sprintf('%s',data_filename);
[s FRdomain] = tfdomain(200,ROI);
% URM_tfunc = URM_freq_response(data_filename,[],[],FRdomain);
[A E B C] = realization(data_filename);
[multH R] = make_SI_op(A,E,B,s0);
N = size(A,1);
tf_err2 = 1;
m = size(R,2);

result = band_Arnoldi(m,R,multH);  % Process the initial band
n = size(result.V,2);

last_ROM_tfunc = [];
j = 0;
while tf_err2 > tfunc_tol
	result = band_Arnoldi(n+stepsize,R,multH,[],n,result);
	
	n = n + stepsize;
	j = j + 1;
	V = make_basis_real(result.V);
	nreal = size(V,2);
	ROM_tfunc = transfer_function(V,A,E,C,B,s);
	
	if ~isempty(last_ROM_tfunc)
		iter(j) = n;
% 		tf_err1(j) = tfunc_err(URM_tfunc,ROM_tfunc);
tf_err1(j) = inf;
		tf_err2(j) = tfunc_err(ROM_tfunc, last_ROM_tfunc);
		fprintf('%d: %g and %g\n',n,tf_err1(j),tf_err2(j));
	end
    last_ROM_tfunc = ROM_tfunc;
end
fprintf('\n');
semilogy(iter,tf_err1','.', iter,tf_err2','+');
legend('rel-error','rel-diff')
title(model_name);

% count flops
	if ~isreal(s0)
		flops = 4 * result.flops;
	else
		flops = result.flops;
	end
	

fprintf('\niterations: %d,  ROM size: %d,  rel-err: %g, flops: %d\n',n,nreal,tf_err1(j),flops)


[VROM VROM_split] = make_basis_real(result.V);
 singVsplit = svd(VROM_split);
 er = rank(VROM_split);
figure;
semilogy(singVsplit,'.');
title(sprintf('%s: svd(V_{ROM}^*), LI: %d',data_filename,er/size(VROM_split,2) ))


	function [err fullerr] = tfunc_err(URM_FR,ROM_FR)
		Ls = size(URM_FR,3);
		fullerr = zeros(1,Ls);
		for i = 1:Ls
			fullerr(i) = norm(URM_FR(:,:,i) - ROM_FR(:,:,i)) / norm(URM_FR(:,:,i));
		end
		err = norm(fullerr);
	end

end
