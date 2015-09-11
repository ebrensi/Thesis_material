function y = MultAr2(x)
%  Computes y = Ar * x
%     implicitly using complex matrix-vector product
%    

global A;

N = length(x)/2;
xc = complex(x(1:N), x(N+1:end));
tmp = A * xc;
y = [real(tmp); imag(tmp)];

