function slvals = lscale(vals) 
% lvals = lscale(vals) 
% 
% scale the elements of an array vals to their logs.  This is to be used 
%  for values of large magnitude, as anything smaller than 1 is rounded to
%  zero.

svals = sign(vals);
lvals = log10(abs(vals));
lvals((lvals < 0)) = 0;
slvals = svals.*lvals;

