function [An En Cn Bn] = project(V,A,E,C,B)
% [An En Cn Bn] = project(V,A,E,C,B)
%
% perform the orthogonal projection of realization (C,G,c,b) on to space spanned by V
An = V'*A*V;
En = V'*E*V;
Bn = V'*B;
Cn = V'*C;

