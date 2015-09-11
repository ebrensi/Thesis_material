% generates plots for Freund's talk(s) at Monterey Linear Algebra conference.  
close all;

%% 1. real s0 close to i-axis:
 AT5('ex308s_1',30,[],(1e-2)*1e10)

%%  2a. s0 moved right on the real-axis
AT5('ex308s_1',30,[],(2.5)*1e10);

%% 2b.  s0 moved up on i-axis 
AT5('ex308s_1',30,[],(1e-2 +2.5i)*1e10)

%% 3. multiple s0 in C produce better model than one in R
AT5('ex1841s_3',80,[],pi*1e10);     % rel_err = 0.0117 
AT5('ex1841s_3',[7 7 7],1e-4,[.8 1 1.5]*pi*1i*1e10); % rel_err = 0.0118

%% 4. single point s0 in C that implies a model half the size of a single
%% point in R.
AT5('ex1841s_3',100,[],pi*1e10);  % n=100 at pi*1e10 -> rel_err = 2.1195e-004
k = 29; s0 = 9.9208e+009 +3.3069e+010i;
AT5('ex1841s_3',k,[],s0);  % n=58 -> rel err = 2e-004


AT5('ex308s_2',100,[],pi*1e10); % rel_err = 1e-4
AT5('ex308s_2',[18 18 18],1e-3,(.3+[1 3.6 6]*2*pi*1i)*1e9);

%% 5. Adaptive method 
[v e] = AT5('ex1a',94,1e-3,pi*1e10);  % n = 94, rel_err = 1.6247e-004

defl_tol = sqrt(eps);  
convergence_tol = 1e-4; 
almost_cnvrgd_tol = convergence_tol; 
wt_good = 0.98;  
s0 = (1+1i)*pi*1e10;
[v e] = AT5a('ex1a',[20 5 5 4]); % n = 57, rel_err = 9.9432e-005
