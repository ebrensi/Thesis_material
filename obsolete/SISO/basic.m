function [Hk Vk] = basic(data,k,s0,flip)
% [Hk Vk] = basic(data,k,s0,flip)
% this is straight-forward full arnoldi to compute the SISO transfer function
%  via explicit projection and via hessian matrix Hk.

if nargin < 4
    flip = false;
    if nargin < 3
        s0 = pi * 1e10;
    end
end

[C G c b] = inputdata(data,flip);

tfunc_real = abs(tfunc_urm(data,C,G,c,b));

[A r0] = makeA(C,G,b,s0);

[Vk Hk] = arnoldi(A,r0,k);
tfunc = tf_hess(Vk,Hk,r0,c,s0);
tfunc2 = tf_proj(Vk,C,G,c,b);

ntf = norm(tfunc_real,inf);
err_h = norm(tfunc_real(:) - tfunc(:),inf)/ntf;
err_p = norm(tfunc_real(:) - tfunc2(:),inf)/ntf;

% plot it on a log scale
% 	h = figure('visible','off');
loglog(f,tfunc_real,'r',f,tfunc,f,tfunc2)
xlabel('f','fontsize',12);
ylabel('|H_k(s)|','fontsize',12,'Rotation',90');
title(sprintf('%s:  k = %d  \t flip: %d',data,k,flip));
legend('actual',sprintf('hess: %3.4g',err_h),sprintf('PRIMA: %3.4g',err_p));
% 	saveas(h,sprintf('A1_%s_%d.png',data,k));
%close h


%% --------------------------------------------------------------

end
