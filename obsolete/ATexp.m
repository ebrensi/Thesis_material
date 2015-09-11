function ATexp(N,m,nvals)
% ATexp(N,m,nvals)
%
% N = size of test matrix to be generated randomly
% m = MIMO (# of input-output variables)
% nvals = # iterations for run(s)


%% Internally defined parameters
defl_tol = sqrt(eps);  % deflation tolerance for band_Arnoldi
convergence_tol = 1e-6; % poles (eigs) with this relative residual are considered converged

%% Setup
rmax = size(nvals,2);

[Amat Armat A Ar R Rr0] = testmat(N);
c = ones(N,1);
cr = [c; -1i*c];
s0 = 0;

%%%% 
tfunc_hess_eqrl = tf_hess([],Armat,Rr0,cr,s0);
tfunc_hess_cplx = tf_hess([],Amat,R,c,s0);
chk = abs([tfunc_hess_cplx; tfunc_hess_eqrl])';
loglog(chk,'.')
legend('cplx','eqrl');
%%%%%%%%%

Nr = 2*N;
n = nvals;
n_tot = sum(nvals);

Vr = zeros(Nr,n_tot*m);
Yr = zeros(Nr,n_tot*m);
Y = [];

yk = 0;
r = 0;
k_tot = 0;
%% Thick Restart band Arnoldi iterations
while r < rmax
    % generate n Arnoldi vectors,
    if r == 0
        %  First time there is no data from previous runs
        tol.defl_flag = 1;
        tol.defl_tol = defl_tol;
        tol.normA_flag = 1;
        tol.normA = 0;
        n0 = 0;
        result = band_Arnoldi(n(1),Rr0,Ar,tol);
        
        if result.exh_flag > 0
            error('Cannot continue! abort after run %d',r);
        end
    else
        clear result0 tol Iv
        % Subsequent runs avoid re-discovering information in
        % the subspace Y, which we already know is invariant.
        Rr0 = [vnext Yr(:,1:yk)];
        tol = result.tol;
        result = band_Arnoldi(n(r+1),Rr0,Ar,tol);
    end
    
    H = result.H;
    vnext = result.Vh_defl;  % !!!possibly remove nearly-zero vectors here!!!
    Vtmp = result.V;
    V = Vtmp(:,n0+1:result.n0);
    %  V is the set of equivalent real vectors obtained from the
    %  previous run of Arnoldi.
    
    vm = size(V,2);
    
    % combine them with the vectors in storage
    Vr(:,k_tot+1:k_tot+vm) = V;
    k_tot = k_tot+vm;
    
    r = r+1;
    if r < rmax
        % Extract a spanning set of Y of (nearly) invariant subspace of V
        Y = deflate(V,H,vnext,convergence_tol);
        l = size(Y,2);
        if l > 0
            % add Y to the (nearly) invariant subspace already discovered
            Yr(:,yk+1:yk+l) = Y;
            yk = yk+l;
        end
    end
end

if k_tot < size(Vr,2)
    Vr(:,k_tot+1:end) = [];
end
if yk < size(Yr,2)
    Yr(:,yk+1:end) = [];
end


V = orth(reshape(Vr,N,[]));


%% -------------------------------------------------------------
% Deflates the set of vectors V, which is provided in equivalent-real
% form. The invariant space Y is returned in equivalent-real form.
    function Y  = deflate(V,H,vnext,convergence_tol)
        [mu W rr tft L] = approx_poles_hess(V,H,vnext,[real(c);imag(c)],Rr0,0);
        
        cvrgd = rr <= convergence_tol;
        n_cvrgd = nnz(cvrgd);
        
        fprintf('%d converged ritz vecs\n',n_cvrgd);
        
        if n_cvrgd > 0
            Wgood = W(:,cvrgd);
            Zgood = V*Wgood;
            rrgood = rr(cvrgd);
            Lgood = L(cvrgd);
            
            % explicitly determine residual error 
            %  this should correspond to rr determined via Arnoldi decomp
            ArZ = Armat*Z;
            ZL = Z*diag(L);
            rr_exp = max(abs(ArZ - ZL))' ./ max(abs(ArZ))';
            rrgood_exp = rr_exp(cvrgd);
            
            chk = [L log(rr) log(rr_exp)];
            
            Zt = Zgood(1:N,:);
            Zb = Zgood(N+1:end,:);
            
            Xe =  (Zt + 1i*Zb ) / sqrt(2);
            actual_eig_idx = sum(abs(Xe).^2)>.5;
            Xe = Xe(:,actual_eig_idx);
            
            Zgood = Zgood(:,actual_eig_idx);
            rrgood = rrgood(actual_eig_idx);
            rrgood_exp = rrgood_exp(actual_eig_idx);
            Lgood = Lgood(actual_eig_idx);
            
            Zalt = [Xe ; -1i*Xe] / sqrt(2);
            
            % Now we compare Zgood with Zalt
            ArZalt = Armat*Zalt;
            rralt = max(abs(ArZalt - Zgood*diag(Lgood)))' ./ max(abs(ArZalt))';
            chk = [rrgood_exp rrgood rralt];
%             semilogy(chk)
%             legend('arnoldi','explicit','reformed')
%             
            Y = [real(Zgood) imag(Zgood)];
        else
            Y = [];
        end
    end

end % of main function
