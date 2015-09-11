function d = dst2s(mu,s)
% compute didtance of a pole to the interval of relevant s values 

smin = min(imag(s));
smax = max(imag(s));

d = zeros(length(mu),1);
below = imag(mu) < smin;
above = imag(mu) > smax;
inrange = ~(above | below); 
d(inrange) = abs(imag(mu(inrange)));

d(below) = abs(mu(below)-smin);
d(above) = abs(mu(above)-smax);
end
