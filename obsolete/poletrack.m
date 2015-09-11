function poletrack(H,vnext,s0,pole)

m = size(H,1);
eta = norm(vnext);

err = inf(m);
rderr = inf(m);
pole = zeros(m);

fprintf('computing poles...\n')
for j = 1:m
    [mu rr] = getpoles(j);
    err(j,1:j) = rr;
    pole(j,1:j) = mu;
    if j > 1
       pmu = pole(j-1,1:j-1).';
       perr = err(j-1,1:j-1).';
       [rd  ix] = mindist(mu,pmu);
       rderr(j,1:j) = rd;
    end
end

err(err<eps) = eps;
rrfix = err(:);
rdfix = rderr(:);
good = rrfix<1e-5;
loglog(rrfix(good),rdfix(good),'.');
xlabel('rr'); ylabel('rd')
    function [mu rr] = getpoles(j)
        [W D] = eig(H(1:j,1:j));
        lambda = ordeig(D);
        rr = eta * abs(W(end,:)).' ./ abs(lambda);
        mu = s0 + 1./lambda;
    end
end

