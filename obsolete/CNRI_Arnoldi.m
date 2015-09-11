function result = CNRI_Arnoldi(nmax,R,MultA,tol,n0,result0)
%
%  Band Arnoldi process
%  (with a quick hack to implement CNRI) 
%
%  This implementation is essentially Algorithm 6.1 in
%  
%  Roland W. Freund, Model reduction methods based on Krylov subspaces, 
%  Acta Numerica, 12 (2003), pp. 267-319.
% 
%  Usage:  result = band_Arnoldi(nmax,R,'MultA')
%          result = band_Arnoldi(nmax,R,'MultA',tol)
%          result = band_Arnoldi(nmax,R,'MultA',tol,n0,result0)
%          
%       nmax = number of Arnoldi vectors to be generated
%          R = matrix the columns of which are the initial vectors
%      MultA = name of a function that performs matrix-vector 
%                     multiplication with the square matrix A
%
%      It is assumed, but not checked, that R and A have the same
%      number of rows.
%
%      If tol is not provided as an input, the following default
%      values are used:
%
%              tol.defl_flag = 1
%              tol.defl_tol = sqrt(eps)    (where eps is machine precision)
%              tol.normA_flag = 1
%
%      If n0 is not provided as an input, the function uses
%                           n0 = 0
%
%  On return:  result is a structure that contains (among othe quantities)
%
%         result.V = matrix V the columns of which are the first nmax 
%                                                  Arnoldi vectors  
%
%         result.H = the nmax x nmax matrix H that represents the projection of 
%                    the matrix A onto the space spanned by the columns of V;
%                    H, V, and A are connected via the relation H = V' * A * V.
%
%       result.rho = the matrix rho that contains the coefficients used to turn
%                    the initial vectors (in R) into the first Arnoldi vectors;
%                    rho, V, and R are connected via the relation rho = V' * R.
% 
%       result.CNRI = a vector of CNRI values for this cycle.
%
%  If the structure tol is provided as an input, it needs to contain 
%
%    tol.defl_tol = unscaled deflation tolerance
%
%  and values of the following to flags:
%
%    tol.defl_flag = 0  use unscaled deflation tolerance
%                    1  use scaled deflation tolerance
%                    2  use scaled deflation tolerance only
%                       for initial block R
%                    3  use scaled deflation tolerance except
%                       for initial block R
%
%    tol.normA_flag = 1  an estimate for the norm of A is
%                        generated within the algorithm
%                     2  an estimate for the norm of A is
%                        provided as tol.normA
%
%  If tol.normA_flag = 2, then an an estimate for the norm of A is
%                      needs to be provided as 
%
%    tol.normA
%
%  If the structure "result0" is provided as an input, then "result0"
%  needs to be the output structure of an earlier call to band_Arnoldi;
%  In this case, the input "tol" is ignored, and the routine uses 
%  tol = result0.tol.
%
% ***********************************************************
%
%  This routine can be run in incremental fashion.
%  Here are two examples.
%
%  Example 1:   result0 = band_Arnoldi(n0,R,'MultA',tol)
%                result = band_Arnoldi(nmax,R,'MultA',tol,n0,result0)
%
%  *  The first call of "band_Arnoldi" runs the band Arnoldi process 
%  *  for n0 step.
%  *
%  *  The second call of "band_Arnoldi" resumes the iteration
%  *  at step n0 + 1 and runs it until step nmax.
%  * 
%  *  Of course, this usage assumes that nmax > n0, and this is 
%  *  checked inside the routine.  
%
%  Example 2 (Arnoldi process, run one step at a time):
%
%      result0 = band_Arnoldi(1,R,'MultA',tol,0,[]);
% 
%      for n = 1 : nmax - 1,
%  
%        result = band_Arnoldi(n+1,R,'MultA',[],n,result0);
%   
%        result0 = result;
%   
%      end
%  
%  *  This will run the Arnoldi process for nmax steps.
%
% ***********************************************************
%
%  Roland W. Freund
% 
%  last change:  July 13, 2005
%
% ***********************************************************
%
if nargin < 3,
  error('** Not enough input arguments! **')
end
%
if nargin == 3,
   tol.defl_flag = 1;
   tol.defl_tol = sqrt(eps);
   tol.normA_flag = 1;
end
%
if nargin <= 4,
  n0 = 0;
end
%
if nmax <= n0,
   error('** nmax is not large enough;  need to have nmax > n0 **')
end  
%
if n0 > 0,
  V = result0.V;
  Vh_defl = result0.Vh_defl;
  rho = result0.rho;
  H = result0.H;
  m = result0.m;
  mc = result0.mc;
  Iv = result0.Iv;
  n0_check = result0.n0;
  tol = result0.tol;
  exh_flag = result0.exh_flag;
%  
  if exh_flag > 0,
    fprintf(' \n')
    disp('**-----------------------------------------------------**')    
    disp('** previous run ended due to exhausted Krylov subspace **')    
    disp('**-----------------------------------------------------**')
    fprintf(' \n') 
    result = result0;    
    return
  end  
%
  if n0 ~= n0_check,
     error('** n0 does not match the value of n0 in result0 **')
  end
%  
  n1 = n0 + 1;
%
else
  [dummy,m] = size(R);
  V = zeros(dummy,0);
  Vh_defl(:,1:m) = R;
  rho = [];
  H = [];
  mc = m;
  Iv.ph = 1:m;
  Iv.I = [];
  Iv.pd = [];  
  Iv.nd = 0;
  n1 = 1;
  exh_flag = 0;  
end
%
result.m = m;
%
%  Tolerances and flags for deflation
%
defl_tol = tol.defl_tol;
defl_flag = tol.defl_flag;
normA_flag = tol.normA_flag;
%
if (normA_flag == 1) && (n0 == 0)
  normA = 0;  % change this.  We might have an estimate for normA despite n0==0.   
else
  normA = tol.normA;
end
%
CNRI = zeros(nmax-n1+1,1);
for n = n1 : nmax,
%
%  Build n-th Arnoldi vector v_n
%
%  Step (1)  (If necessary, deflate v vector)
%
  foundvn = 0;
%
  while foundvn == 0,
%
    [mc,foundvn,Vh_defl,Iv,normv] = deflation(0,n,m,mc,foundvn, ...
                               Vh_defl,R,Iv,defl_flag,defl_tol,normA);
%
    if mc == 0,
      disp('**---------------------------------------**')
      disp('** There are no more Krylov vectors, and **')
      disp('** so the process has to terminate: STOP **')
      disp('**---------------------------------------**')
      disp(['  Number of Arnoldi steps performed: ' num2str(n-1)])
      result.V = V;
      result.Vh_defl = Vh_defl;
      result.rho = rho;
      result.H = H;
      result.mc = mc;
      result.Iv = Iv;
      result.Vh_defl = Vh_defl;
      tol.normA = normA;
      result.tol = tol;
      result.n0 = n - 1;
      result.exh_flag = 1;
	  result.CNRI = CNRI;
      return
    end
%        
%    End of: while foundvn == 0
%
  end
%
%   Make sure rho has n rows
%
  rho(n,1) = 0;
%
%   Step (2)  (Normalize v_n)
%
  V(:,n) = Vh_defl(:,Iv.ph(1)) / normv;
  if n > mc,
    H(n,n-mc) = normv;
  else
    rho(n,n-mc+m) = normv;
  end
%
%   Step (3)  (Orthogonalize the candidate vectors against v_n)
%
  ivph1 = Iv.ph(1);
  for k = 1 : mc -1,
    ikp1 = Iv.ph(k+1);
    tmp = (V(:,n))' * Vh_defl(:,ikp1);
    Vh_defl(:,ikp1) = Vh_defl(:,ikp1) - V(:,n) * tmp;   
    if k > mc - n,
      H(n,k-mc+n) = tmp;
    else
      rho(n,k-mc+n+m) = tmp;
    end
    Iv.ph(k) = ikp1;
  end
  Iv.ph(mc) = ivph1;
%
%   Step (4)   (Advance block Krylov subspace)
%
%   Step (4a)  (Compute tmpv = A * V(:,n))
%
  tmpv = MultA(V(:,n));
  V1 = [V tmpv];
  vals = svd(V1);
% CNRI(n) = vals(n+1)/vals(1);
  CNRI(n) = 1/cond(V1,2);
%
  if normA_flag == 1,
    normA  = max( [ normA, norm(tmpv,2) ] );
  end
%
%   Step (4b)  (Orthogonalize tmpv against previous Arnoldi vectors)
%
  for k = 1 : n,
    H(k,n) = (V(:,k))' * tmpv;
    tmpv   = tmpv - V(:,k) * H(k,n);
  end
%
  Vh_defl(:,Iv.ph(mc)) = tmpv;
%
%   Step (5)  (Compute entries in H resp. rho due to nonzero deflated vectors)
%    
  nd = Iv.nd;
%    
  for k = 1 : nd,
    tmp = (V(:,n))' * Vh_defl(:,Iv.pd(k));
    ivk = Iv.I(k);
    if ivk > 0,
      H(n,ivk) = tmp;
    else
      rho(n,ivk + m) = tmp;
    end 
  end
%
  result.n0 = n;
%
end
%
result.V = V;
result.Vh_defl = Vh_defl;
result.rho = rho;
result.H = H;
result.mc = mc;
result.Iv = Iv;
tol.normA = normA;
result.tol = tol;
result.exh_flag = 0;
result.CNRI = CNRI;
%
