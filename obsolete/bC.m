function [Hk Vk] = bC(data,k,s0)
% [Hk Vk] = b2(data,k,s0)
% this is straight-forward full complex arnoldi to compute the SISO
% transfer function via explicit projection.  Uses Freund's band_arnoldi code.

if nargin < 3
    s0 = (1+1i)* pi * 1e10;
end

[C G c b] = inputdata(data);
[N m] = size(b);
tfunc_unreduced = abs(tfunc_urm(data,C,G,c,b));

[MultA R0] = makeA(C,G,b,s0);

tic;
result = band_Arnoldi(k,R0,MultA);
Vc = result.V;
% Vs = reshape([real(Vc); imag(Vc)],N,[]);
% Vk = orth(Vs);
Vk = Vc;
toc 

tfunc = tf_proj(Vk,C,G,c,b);

rel_err = plot_tfunc(tfunc,tfunc_unreduced,'reduced','unreduced');
rel_err

Ck = Vk'*C*Vk;
Gk = Vk'*G*Vk;
if ~isposreal2(Ck,Gk)
    fprintf('reduced model not pos real\n');
else
    fprintf('reduced model pos real\n');
end
% [mu Z rr tf_terms] = approx_poles_proj(Vk,C,G,c,b,s0,false);
%  wt = tf_terms_weight(tf_terms);
% plot_poles(mu, wt,rr,s0);

end
