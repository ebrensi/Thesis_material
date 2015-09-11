function result = paired_band_Arnoldi_mod(nmax,R,MultAr,tol,n0,result0)
%
%  Paired Band Arnoldi process
%  
%  Usage:  result = band_Arnoldi(nmax,R,'MultAr')
%          result = band_Arnoldi(nmax,R,'MultAr',tol)
%          result = band_Arnoldi(nmax,R,'MultAr',tol,n0,result0)
%          
%       nmax = number of Arnoldi vectors to be generated
%          R = matrix the columns of which are the initial vectors
%     MultAr = name of a function that performs matrix-vector 
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
%              tol.defl_tol_gs = sqrt(eps);
%              tol.normA_flag = 1
%
%      If n0 is not provided as an input, the function uses
%                           n0 = 0
%
%  On return:  result is a structure that contains (among othe quantities)
%
%         result.V = matrix V the columns of which are the first nmax 
%                                                     Arnoldi vectors  
%         
%         result.W = matrix W the columns of which are a real
%                    orthonormal basis of the paired Krylov subspace
%
%         result.H = the nmax x nmax matrix H that represents the projection of 
%                    the matrix A onto the space spanned by the columns of V;
%                    H, V, and A are connected via the relation H = V' * A * V.
%
%       result.rho = the matrix rho that contains the coefficients used to turn
%                    the initial vectors (in R) into the first Arnoldi vectors;
%                    rho, V, and R are connected via the relation rho = V' * R.
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
%  Example 1:   result0 = band_Arnoldi(n0,R,'MultAr',tol)
%                result = band_Arnoldi(nmax,R,'MultAr',tol,n0,result0)
%
%  *  The first call of "paired band_Arnoldi" runs the band Arnoldi process 
%  *  for n0 step.
%  *
%  *  The second call of "paired_band_Arnoldi" resumes the iteration
%  *  at step n0 + 1 and runs it until step nmax.
%  * 
%  *  Of course, this usage assumes that nmax > n0, and this is 
%  *  checked inside the routine.  
%
%  Example 2 (Arnoldi process, run one step at a time):
%
%      result0 = paired_band_Arnoldi(1,R,'MultAr',tol,0,[]);
% 
%      for n = 1 : nmax - 1,
%  
%        result = paired_band_Arnoldi(n+1,R,'MultAr',[],n,result0);
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
%  last change:  May 11, 2010
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
   tol.defl_tol_gs = sqrt(eps);
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
  W = result0.W; 
  rho = result0.rho;
  H = result0.H;
  m = result0.m;
  mc = result0.mc;
  lc = result0.lc;
  Iv = result0.Iv;
  Iw = result0.Iw;
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
  W = zeros(dummy/2,0);
  Vh_defl(:,1:m) = R;
  Vcount = [];
  Wcount = [];
  rho = [];
  H = [];
  mc = m;
  lc = 0;
  Iv.ph = 1:m;
  Iv.I = [];
  Iv.pd = [];  
  Iv.nd = 0;
  n1 = 1;
  Iw.I = [];
  Iw.nd = 0;
  exh_flag = 0;  
end
%
result.m = m;
%
%  Tolerances and flags for deflation
%
defl_tol = tol.defl_tol;
defl_tol_gs = tol.defl_tol_gs;
defl_flag = tol.defl_flag;
normA_flag = tol.normA_flag;
%
if (normA_flag == 1) && (n0 == 0),
  normA = 0;
else
  normA = tol.normA;
end
%
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
      result.W = W;
      result.Vcount = Vcount;
      result.Wcount = Wcount;
      result.Vh_defl = Vh_defl;
      result.rho = rho;
      result.H = H;
      result.mc = mc;
      result.lc = lc;
      result.Iv = Iv;
      result.Iw = Iw;
      tol.normA = normA;
      result.tol = tol;
      result.n0 = n - 1;
      result.exh_flag = 1;
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
%  Orthogonalize the top half of V(:,n) against the 
%  previous w_j's

   N = length(V(:,n))/2;
   tmpw = V(1:N,n);
   for j = 1:lc,
     tmp = W(:,j)' * tmpw;
     tmpw = tmpw - W(:,j) * tmp;
   end
%
   norm_tmpw = norm(tmpw);
   if norm_tmpw > defl_tol_gs,
     lc = lc + 1;
%      size(tmpw)
     W(:,lc) = tmpw / norm_tmpw;
   else
     Iw.nd = Iw.nd + 1;
     Iw.I = [Iw.I n];
   end
%
%  Orthogonalize the bottom half of V(:,n) against the 
%  previous w_j's
%
   tmpw = V(N+1:end,n);
   for j = 1:lc,
     tmp = W(:,j)' * tmpw;
     tmpw = tmpw - W(:,j) * tmp;
   end
%
   norm_tmpw = norm(tmpw);
   if norm_tmpw > defl_tol_gs,
     lc = lc + 1;
     W(:,lc) = tmpw / norm_tmpw;
   else
     Iw.nd = Iw.nd + 1;
     Iw.I = [Iw.I -n];    
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
  tmpv = feval(MultAr,V(:,n));
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
  Vcount(n) = size(V,2);
  Wcount(n) = size(W,2);
%
end
%
result.V = V;
result.W = W;
result.Vcount = Vcount;
result.Wcount = Wcount;
result.Vh_defl = Vh_defl;
result.rho = rho;
result.H = H;
result.mc = mc;
result.lc = lc;
result.Iv = Iv;
result.Iw = Iw;
tol.normA = normA;
result.tol = tol;
result.exh_flag = 0;
%
