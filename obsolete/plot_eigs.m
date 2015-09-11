function plot_eigs(ritz_vals, weights,rr,axs)
% -- Plot eigenvalues and ritzvals--
maxsize = 96;

if exist('A_eigs','var');
    scatter(real(A_eigs),imag(A_eigs),'k.');
end

% hold on

if nargin == 1
    scatter(real(ritz_vals),imag(ritz_vals),10,'r','filled');
elseif nargin == 2
    sz = weights/max(weights) * maxsize;
    scatter(real(ritz_vals),imag(ritz_vals),sz,'r','filled');
elseif nargin == 3
%     colormap('bone');
    sz = weights/max(weights) * maxsize;
    scatter(real(ritz_vals),imag(ritz_vals),sz,-log(rr),'filled');
end
colorbar;

if exist('axs','var')
    axis(axs);
end
hold off

end
