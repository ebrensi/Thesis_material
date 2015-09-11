%  Example driver for paired_band_Arnoldi.m
%    
%  This driver creates a random N x N matrix A
%  and a random matrix R of m = 5 initial vectors;
%  R is chosen such that some deflations will occur
%  early on in the band Arnoldi process.
%
%  It then runs the band Arnoldi process for nmax steps.
%  
%  * The first run is in incremental mode, with one Arnoldi step
%    run at a time.
%
%  * In the second run, the band Arnoldi process is run for all
%    nmax steps via a single call.
%
% ***********************************************************
%
%  Roland W. Freund
%
%  last change:  May 11, 2010
%
% ***********************************************************
%
global A;

tol.defl_tol = sqrt(eps);
tol.defl_flag = 1;
tol.normA_flag = 1;
tol.defl_tol_gs = sqrt(eps);

m = 5;

N = input('Size N of the N x N matrix A = ')

nmax = input('Number nmax of Arnoldi vectors to be generated = ')

rand('seed',0);
A         = rand(N,N);
A         = A + 1i * rand(N,N);
Rc        = rand(N,m) - 0.5 * ones(N,m);
Rc        = Rc + 1i * (rand(N,m) - 0.5 * ones(N,m));
Rc(:,2)   = Rc(:,1) + 1e-12 * rand(N,1);
Rc(:,m)   = A*A*A*Rc(:,1) + 1e-12 * rand(N,1);
Rc(:,m-1) = A*A*Rc(:,1) + 1e-12 * rand(N,1);

R = [real(Rc); imag(Rc)];

fprintf(' \n')
disp('**-------------------------------**')
disp('** First run in incremental mode **')
disp('**-------------------------------**'), fprintf(' \n')

result0 = paired_band_Arnoldi(1,R,'MultAr',tol,0,[]);

for n = 1 : nmax,

  result_i = paired_band_Arnoldi(n+1,R,'MultAr',[],n,result0);
  
  result0 = result_i;
  
end  

fprintf(' \n')
disp('**------------------------------------**')
disp('** Second run: all nmax steps at once **')
disp('**------------------------------------**'), fprintf(' \n')

result = paired_band_Arnoldi(nmax+1,R,'MultAr',tol,0,[]);

