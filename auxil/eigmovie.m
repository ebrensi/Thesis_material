function mov = eigmovie(H,vnext,s0,convergence_tol,ell,ax)
% mov = eigmovie(mov = eigmovie(H,vnext,s0,convergence_tol,ell,ax,mov))
% 
% Returns an animated visualization of eigenvalues of leading principle 
%  submatrices of H

m = size(H,1);
eta = norm(vnext);

[mu_m rr_m] = getpoles(m);

if ~exist('ax','var') || isempty(ax)
    % Start with the largest matrix to determine axis size
    max_i = max(imag(mu_m));
    min_i = min(imag(mu_m));
    max_r = lscale(max(real(mu_m)));
    min_r = lscale(min(real(mu_m)));
    ax = [min_r max_r min_i max_i];
end

if nargout
    movie_index = 1;
end


set(gca,'nextplot','replacechildren');

for j = ell+1:m
    if j < m
        [mu rr] = getpoles(j);
    else
        mu = mu_m;
        rr = rr_m;
    end
    
    
    % plot the poles from this run
    plot_poles(mu,[],rr,s0,true,convergence_tol);
    axis(ax);
    title(sprintf('dim(Y)=%d,  j=%d',ell,j-ell));
    drawnow
    if nargout
        mov(movie_index) = getframe(gcf);
        movie_index = movie_index+1;
    end
end


    function [mu rr] = getpoles(j)
        % compute approximate poles as eigenvalues of Arnoldi matrix
        [W D] = eig(H(1:j,1:j));
        lambda = ordeig(D);
        rr = eta * abs(W(end,:)).' ./ abs(lambda);
        mu = s0 + 1./lambda;
    end
end % of main function 
