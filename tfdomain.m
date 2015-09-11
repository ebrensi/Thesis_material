function [s frq] = tfdomain(npts,freq_interval,opt)
% [s frq] = domain_tfunc(npts,freq_interval)
%
% frq and s are using in so many routines and they never change, so 
%  I wrote a function to load them in memory.  If we ever change the 
%  desired frequency range, we can then do it globally via modifying this
%  function.
%  
%  Default is 3200 points from 10^9 to 10^10, distributed logly



if ~exist('npts','var') || isempty(npts)
	npts = 3200;                 % default number of sample points
end
if ~exist('freq_interval','var') || isempty(freq_interval)
	freq_interval = [9 10];     % default frequency range 
end
if exist('opt','var') && strcmp(opt,'linear') 
	frq = linspace(10^freq_interval(1),10^freq_interval(2),npts);
else
	frq = logspace(freq_interval(1),freq_interval(2),npts);
end
omega = 2*pi*frq;
s = 1i * omega;

