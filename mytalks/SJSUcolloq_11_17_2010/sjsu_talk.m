% part I
H = [1 2 3; 4 5 6; 2 7 3]
r = [1 1 1].'

% part II
S = [ 1 2 3; 1 2 2; 1 3 1]
M = diag([2 .5 1])
A = S*M*inv(S)
 
% r = [1 1 1; 3 -6 0]'

r = [7 -5 11]';
for j = 1:10
    r = A*r;
    j
    latexmatrix2(sum(r,2),'%5.4g')
end

 