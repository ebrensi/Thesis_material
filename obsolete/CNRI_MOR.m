function [Hn Vn] = CNRI_MOR(data,nmax,s0,flip)
% [Hn Vn] = basic(data,nmax,s0,flip)
% this is straight-forward full arnoldi to compute the SISO transfer function
%  via explicit projection and via hessenberg matrix Hk.

freund = 0;

if nargin < 4
    flip = false;
    if nargin < 3
        s0 = pi * 1e10;
    end
end

[C G c b] = inputdata(data,flip);

URM = abs(tfunc_urm(data,C,G,c,b));

[A r0] = makeA(C,G,b,s0);

if freund == 1
	result = CNRI_Arnoldi(nmax,r0,A);
	Vn = result.V;
	Hn = result.H;
	rho = result.rho;
	CNRI = result.CNRI;
else
	[Vn v_next Hn h_next] = arnoldi(A,r0,nmax);
end
Hn_subdiag = diag(Hn,-1);
n = size(Vn,2);
rho = zeros(n,1);
rho(1) = norm(r0);

MOR_H = abs(tf_hess(Vn,Hn,rho,c,s0));
MOR_P = abs(tf_proj(Vn,C,G,c,b));

ntf = norm(URM,inf);
err_h = norm(URM(:) - MOR_H(:),inf)/ntf;
err_p = norm(URM(:) - MOR_P(:),inf)/ntf;


% plot it on a log scale
% 	h = figure('visible','off');
[s f] = getS();
loglog(f,URM,'r',f,MOR_H,f,MOR_P)
xlabel('f','fontsize',12);
ylabel('|H_k(s)|','fontsize',12,'Rotation',90');
title(sprintf('%s:  n = %d  \t flip: %d',data,n,flip));
legend('actual',sprintf('H: %3.4g',err_h),sprintf('proj: %3.4g',err_p));
% 	saveas(h,sprintf('A1_%s_%d.png',data,k));
%close h


if false
	figure
	semilogy(CNRI,'.');
	disp([CNRI [Hn_subdiag;0]])
end

%% --------------------------------------------------------------

end
