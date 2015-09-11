function [Hk Vk] = bR(data,n,nruns,s0)
% [Hk Vk] = bR(data,k,s0)
% this is full eq real arnoldi to compute the SISO
% transfer function via explicit projection.  Uses Freund's band_arnoldi code.

if ~exist('s0','var')
    s0 = pi * 1e10;
end

[C G c b] = inputdata(data);
[N m] = size(b);

tfunc_unreduced = abs(tfunc_urm(data,C,G,c,b));
[Ar Rr0] = makeAr(C,G,b,s0);
Nr = 2*N;
Vr = zeros(Nr,n*m*nruns);

tic;
    if r == 1
        result = band_Arnoldi(n,Rr0,Ar);
        tol = result.tol;
    else
%         result0.V = zeros(Nr,0);
%         result0.Vh_defl(:,1:m) = Rr0;
%         result0.rho = [];
%         result0.H = [];
%         result0.mc = m;
%         Iv.ph = 1:m;
%         Iv.I = [];
%         Iv.pd = [];
%         Iv.nd = 0;
%         result0.Iv = Iv;
%         result0.exh_flag = 0;
%         result0.tol = tol;
        result0 = result;
        
        n0 = size(result.V,2);
        result = band_Arnoldi(n*r,Rr0,Ar,[],n0,result0);
    end
    
end

Vs = reshape(result.V,N,[]);
Vk = orth(Vs);
toc

tfunc = tf_proj(Vk,C,G,c,b);
rel_err = plot_tfunc(tfunc,tfunc_unreduced,'reduced','unreduced');
rel_err

Ck = Vk'*C*Vk;
Gk = Vk'*G*Vk;
if ~isposreal2(Ck,Gk)
    fprintf('ROM is not pos real\n');
else
    fprintf('ROM is pos real\n')
end
%  wt = tf_terms_weight(tf_terms);
% plot_poles(mu, wt,rr,s0);

end
