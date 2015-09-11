function [A E B C frq tfunc] = realization(data)
% [A E B C frq tfunc] = realization(data)

S = load(data);

A = S.A;  
B = S.B;

if isfield(S,'E')
	E = S.E;
else
	E = eye(size(A));
end

if isfield(S,'C')
	C = S.C;
else 
	C = S.B;
end

if size(C,1) ~= size(A,1)
	C = C';
end

if isfield(S,'w') 
	frq = S.w;
else 
	frq = [];
end

if isfield(S,'mag') 
	tfunc = S.mag;
else
	tfunc = [];
end
