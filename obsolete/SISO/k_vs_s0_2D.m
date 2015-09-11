function [s0opt_k s0opt_sz] = opt_s0(datafilename,nkpts,ns0pts)
% [s0opt_k s0opt_sz] = opt_s0(datafilename,nkpts,ns0pts)
%
% determines optimal s0 for (non-restarted) Arnoldi ROM construction
%  within user specified constraints:

err_tol = 1e-2;       % ROM must have relative error less than err_tol
kmax = 130;           % ROM must be computable in kmax iterations or less 
% s0 is located on a specified grid 
rl_rng = linspace(1e-6,2,ns0pts);  % the range of s0 points on the real axis
im_rng = linspace(0,2,ns0pts);     %  " "   on the imag-axis 

%% nkpts evenly spaced k values from 2 to kmax  
kvals = fix(linspace(2,kmax,nkpts));

%% make grid of test expansion points s0 
[X Y] = meshgrid(rl_rng,im_rng);
S0 = (X + Y*1i) *pi * 1e10;  % S0 is the matrix of s0 values

%% Setup 

% import model data
[C G c b] = inputdata(datafilename);
[N p] = size(b);

% get transfer function for the unreduced model to compare against
tfunc_actual = abs(tfunc_urm(datafilename,C,G,c,b))';
ntf = norm(tfunc_actual,inf);

nel = length(S0(:));

% preallocate 
K = zeros(size(S0));
SZ = zeros(size(S0));
ERR = zeros(size(S0));

%% Evaluate ROMs
tic; % used for keeping track of elapsed time
for j = 1:nel
    s0 = S0(j);
    [K(j) SZ(j) ERR(j)] = find_good_model(s0,kvals,err_tol);
    
    et = toc;
    rt = et/j;
    eta = rt*(nel-j);  % rough estimate of how much time is left
    
    if ~isinf(NS(j))
        fprintf('k=%d,\t n=%d\t s0 = (%g + %gi)*pi*1e10\r',K(j),SZ(j),X(j),Y(j));
    end
    fprintf('eta: %d:%02d\r',fix(eta/60),fix(rem(eta,60)));
end
fprintf('\n');

%% Determine best ROMs and plot results 

% Evaluate based on number of iterations (k) of Arnoldi required
%  this gives an idea of how efficiently we can produce the ROM
[min_k im] = min(K(:));
s0opt_k = S0(im);
mesh(X,Y,K);
hold on
plot3(X(im),Y(im),K(im),'r+');
hold off
title(sprintf('%s: opt s0 = (%3.2g, %3.2gi)\t  k = %d,\t sz=%d,\t err = %2.1e',datafilename,X(im),Y(im),K(im),SZ(im),ERR(im)));
xlabel('Re'); ylabel('Im');
% saveas(gcf,sprintf('%s_k',datafilename));  % save the plot as a .fig file


% Evaluate based on order of the ROM
%  this gives an idea of how small we can get our model
[min_sz im] = min(SZ(:));
s0opt_sz = S0(im);
mesh(X,Y,SZ);
hold on
plot3(X(im),Y(im),SZ(im),'r+');
hold off
title(sprintf('%s: opt s0 = (%3.2g, %3.2gi)\t  k = %d,\t sz=%d,\t err = %2.1e',datafilename,X(im),Y(im),K(im),SZ(im),ERR(im)));
xlabel('Re'); ylabel('Im');
% saveas(gcf,sprintf('%s_sz',datafilename));  % save the plot as a .fig file

%% ----------------------------------------------------- 

    function [k sz err]= find_good_model(s0,kvals,err_tol)
   %  Given s0 this function constructs a ROM using k iterations of Arnoldi,
   %    for each k in kvals until it finds a ROM with error less than
   %    err_tol or runs out of kvals. 
   %  If it finds such a 'good' ROM it returns: 
   %  k  = the # of iterations k required to construct the model
   %  sz =  the order of the ROM
   %  err = the error of the ROM
   %
   %  Otherwise everything returned as inf
   
        kvals = sort(kvals);
        [Ar Rr0] = makeAr(C,G,b,s0);
        
        found_good_model = false;
        i = 1;
        while i <= length(kvals) && ~found_good_model
            k = kvals(i);
            if i == 1
                result = band_Arnoldi(k,Rr0,Ar);
            else
                result = band_Arnoldi(k,Rr0,Ar,[],kvals(i-1),result);
            end
            
            Vk = orth(reshape(result.V,N,[]));
            sz = size(Vk,2);
            tfunc_k = abs(tf_proj(Vk,C,G,c,b));
            err = norm(tfunc_actual(:) - tfunc_k(:),inf)/ntf;
            
            if err <= err_tol
                found_good_model = true;
%                 plot_tfunc(tfunc_k,tfunc_actual,sprintf('k=%d, sz=%d',k,sz),'exact');
%                 drawnow; 
            end
        end
        
        if ~found_good_model
            k = inf;
            sz = inf;
            err = inf;
        end
    end
       
end  % end main function
