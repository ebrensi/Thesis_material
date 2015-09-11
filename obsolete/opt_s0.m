function [s0_opt opt_err] = opt_s0(data,kmax,s0_real_init)
% this is an attempt to find the optimum value of s0 via 
%  matlab's built-in optimization routines.

[C G c b] = inputdata(data);
frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;
tfunc_real = abs(tfunc_urm(data,C,G,c,b)');
cplx = [1 sqrt(-1)].'*pi*1e10;

% s0_real_init = [1 1];%[1 1];
[s0_opt_real, opt_err, exitflag,output] = fminsearch(@tfunc_err,s0_real_init);
s0_opt_real
s0_opt = s0_opt_real * cplx;
fprintf('exitflag = %d, s0_{opt} = %g + %gi, err_{opt} = %g\n',exitflag,real(s0_opt),imag(s0_opt),opt_err);

% ----------------------------------------------------------------
    function err = tfunc_err(s0_real)
        s0_real
        s0 = s0_real(:)'*cplx;
        [A r0] = makeA(C,G,b,s0);
        [H V] = arnoldi(A,r0,kmax);
        tfunc = abs(tf_hess(V,H,r0,c,s,s0));
        %     tfunc = tf_proj(V,C,G,c,b,s);
        err = norm(tfunc-tfunc_real,inf)/norm(tfunc_real,inf);
    end
end
