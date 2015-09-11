function H = bode3d

ax = [-6 4 -5 5];

close all
mu = [inf -.2+2i -.2-2i -1 -4+3i]; 
r =  [1, 1, 1, 3,20];

wt = abs(r ./ real(mu))
mx = 20;

[sig omeg] = meshgrid(ax(1):.1:ax(2),ax(3):.1:ax(4));

Y = U(complex(sig,omeg));
Y(abs(Y)>mx) = mx;

% figure;
% mesh(sig,omeg,angle(Y));
% zlabel('arg(H(s))','fontsize',12,'fontweight','demi')
% xlabel('\sigma','fontsize',12,'fontweight','demi')
% ylabel('j\omega','fontsize',12,'fontweight','demi')

figure;

contour(sig,omeg,abs(Y))
xlabel('\sigma','fontsize',12,'fontweight','demi')
ylabel('j\omega','fontsize',12,'fontweight','demi')


figure
mesh(sig,omeg,abs(Y))
zlabel('|H(s)|','fontsize',12,'fontweight','demi')
xlabel('\sigma','fontsize',12,'fontweight','demi')
ylabel('j\omega','fontsize',12,'fontweight','demi')

figure;
omeg = ax(3):.1:ax(4);
H = abs(U(1i*omeg));
plot(omeg,abs(H));
axis([ax(3) ax(4) 0 12]);
ylabel('|H(j\omega)|','fontsize',12,'fontweight','demi')
xlabel('j\omega','fontsize',12,'fontweight','demi')

% figure
% plot(omeg,angle(U(1i*omeg)));
% ylabel('arg(H(j\omega))','fontsize',12,'fontweight','demi')
% xlabel('j\omega','fontsize',12,'fontweight','demi')

    function u = U(s)
        u = 0;
        for i = 1:length(mu)
            if isinf(mu(i))
                u = u+r(i);
            else
                u = u + r(i)./(s - mu(i));
            end
        end
    end
end
