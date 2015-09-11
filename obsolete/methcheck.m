function methcheck(data,kmax)
% this routine compares different methods to implement the matrix-vector 
%   multiplication Av used in the Arnoldi process.

[C G c b] = inputdata(data);
tfunc_actual = abs(tfunc_urm(data,C,G,c,b));

frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
s = 2*pi*sqrt(-1)*frq;
s0 = pi*1e10;

C_s0G = C-s0*G;
r0 = -C_s0G\b;

% method 1 explicit solve
tic;
[H1 V1] = arnoldi(@multA1,r0,kmax);
tfunc1 = tf_hess(V1,H1,r0,c,s,s0);
t1 = toc;
tfunc1_err = tfunc_err(tfunc1);


% method 2 LU p
tic;
[L U p] = lu(C_s0G,'vector');
[H2 V2] = arnoldi(@multA2,r0,kmax);
t2 = toc;
tfunc2 = tf_hess(V2,H2,r0,c,s,s0);
tfunc2_err = tfunc_err(tfunc2);


% method 3 LU PQ
tic;
[L U P Q] = lu(C_s0G);
[H3 V3] = arnoldi(@multA3,r0,kmax);
t3 = toc;
tfunc3 = tf_hess(V3,H3,r0,c,s,s0);
tfunc3_err = tfunc_err(tfunc3);

loglog(frq,tfunc_actual,'k-.',frq,[tfunc1; tfunc2; tfunc3]');
legend('actual','mldivide','LUp','LUPQ');
fprintf('mldivide (%3.2g):\t%g \nLUp (%3.2g):\t\t%g\nLUPQ (%3.2g):\t\t%g\n\n',t1,tfunc1_err,t2,tfunc2_err,t3,tfunc3_err);

    function q = multA1(v)
        q = C_s0G\(G*v);
    end

    function q = multA2(v)
        y = G*v;
        q = U\(L\(y(p)));
    end

    function q = multA3(v)
        y = G*v;
        q = Q * (U\(L\(P*y)));
    end

    function err = tfunc_err(tfunc)
        err = norm(abs(tfunc_actual(:)) - abs(tfunc(:)),inf)/norm(tfunc_actual);
    end

end
