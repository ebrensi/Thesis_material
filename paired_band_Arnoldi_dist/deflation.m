function  [mc,fv,Vh_d,I,nv] = deflation(vw_f,n,m,mc,fv,Vh_d,R,I,d_f,d_t,nA)
%
%  This function checks for deflation, and if necessary, performs
%  the deflation.
%
% ***********************************************************
%
%  Roland W. Freund
%
%  last change:  July 13, 2005
%
% ***********************************************************
%
itmp = I.ph(1);
nv = norm(Vh_d(:,itmp));
%
if n <= mc,
  if (d_f == 1) || (d_f == 2),
    def_fac = norm(R(:,n+m-mc),2);
  else
    def_fac = 1;
  end
else
  if (d_f == 1) || (d_f == 3),
    def_fac = nA;
  else
    def_fac = 1;   
  end
end
%        
act_defl_tol = def_fac * d_t;
%
if ( nv > act_defl_tol )
  fv = 1;
else
%   fprintf(' \n')
%   disp('**--------------------------------------**')
% %
%   switch vw_f
%     case 1
%       disp('** Deflation of next candidate v vector **')
%     case -1
%       disp('** Deflation of next candidate w vector **')
%     case 0
%       disp('** Deflation of next candidate vector **')      
%   end
%  
%   disp('**--------------------------------------**')
%   disp(['  Iteration index (n)    : ' num2str(n)])
%   disp(['  Deflation Tolerance    : ' num2str(act_defl_tol)])
%   disp(['  Norm of deflated vector: ' num2str(nv)])
%

fprintf('v_defl(%d)/H_est  = %g < %g, \t mc = %d\n',n, nv/def_fac,d_t,mc-1);  % 1/15/2014 E. Rensi

%   switch vw_f
%     case 1   
%       disp(['  New right block size   : ' num2str(mc-1)])
%     case -1
%       disp(['  New left block size    : ' num2str(mc-1)])
%    case 0      
%       disp(['  New block size         : ' num2str(mc-1)])    
%   end
%  
%   disp('**--------------------------------------**'), fprintf(' \n')
%
%     Shift the pointers to the candidate vectors
%
  for k = 1 : mc - 1,
    I.ph(k) = I.ph(k+1);
  end
%
%     Update I.nd, I.I, and I.pd
%
  if nv > 0,    
    I.nd = I.nd + 1;
    I.I = [I.I, n - mc];
    I.pd = [I.pd, itmp];
  else
    Vh_d(:,itmp) = zeros(size(Vh_d(:,itmp)));
  end
%
%     Update current block size and resize I.ph
%
  mc = mc - 1;
  I.ph = I.ph(1:mc);
%
end
%        
%    End of: if fv == 0
%


