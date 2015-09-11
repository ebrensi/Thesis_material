function pr = test1(data,s0,alpha,beta)
frq = logspace(8,10,200);

% example1 from later paper
% A = [-2 -1;1 0]; E = eye(2); b = [1 0]';  c = [1 3]';

% % ex1 from original paper
% A = [-3 -2;1 0]; E = eye(2); b = [1 0]'; c = [1 1.5]'; 
% 
% % ex3 from original paper
% A = [0 0 0; 1 0 -4; 0 1 -5]; E = eye(3); b = [6 5 1]'; c = [0 0 1]';  

[A E c b] = inputdata(data);
if ~exist('alpha','var')
    alpha = 0;
end
if ~exist('beta','var')
    beta = 0;
end

% s0 = 100;
n = length(b);
I = speye(n);
A = A - alpha*I;
E = E - beta*I;
pr = isposreal(A,E,c,b,s0);

e1 = eye(length(b),1);


% ------------------------------------------------------------------------

    function H = tfunc(A,E,c,b)
        frq = logspace(8,10,200);
        s = 2*pi*sqrt(-1)*frq;
        Ls = length(s);
        [UA,UE,Q,Z] = qz(full(A),full(E));
        cc = Z'*c;
        bb = Q*b;
        H = zeros(1,Ls);
        optU.UT = true; 
        for j = 1:length(s)
            H(j) = abs(cc'*linsolve(s(j)*UE-UA,bb,optU));
        end
    end % tfunc

end % main function
