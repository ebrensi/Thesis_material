function d = ftest(N,nmax,meth)

m = 1;
A = (rand(N) - 0.5) + 1i* (rand(N) - 0.5);
R = (rand(N,m) - 0.5) + 1i* (rand(N,m) - 0.5);

%% perform complex Arnoldi
% [Vc,Hc] = my_barnoldi(@MultA,R,nmax);
result = band_Arnoldi(nmax,R,@MultA);
Vc = result.V;
Hc = result.H; 

Vc = Vc(:,1:nmax);
Hc = Hc(1:nmax,1:nmax);

Vc_ortho_diff = norm(eye(nmax) - Vc'*Vc,inf);

%% perform equivalent real Arnoldi
Rr = [real(R); imag(R)];
if meth == 1
    % using explict eq real formulation of A
    A1 = real(A);
    A2 = imag(A);
    Ar = [A1 -A2 ; A2  A1];
    
%     [Vr,Hr] = my_barnoldi(@MultAr,Rr,nmax);
    result = band_Arnoldi(nmax,Rr,@MultAr);
    Vr = result.V;
    Hr = result.H;
else
    % using complex matrix-vector product
    result = band_Arnoldi(nmax,Rr,@MultAr2);
    Vr = result.V;
    Hr = result.H;
%     [Vr,Hr] = my_barnoldi(@MultAr2,Rr,nmax);
end

Vr = Vr(:,1:nmax);
Hr = Hr(1:nmax,1:nmax);
Vr_ortho_diff = norm(eye(nmax) - Vr'*Vr,inf);

%%  Split the bases into split sets
Vc1 = zeros(N,0);
Vr1 = zeros(N,0);

for j = 1:nmax,
    Vc1 = [Vc1 real(Vc(:,j)) imag(Vc(:,j))];
    Vr1 = [Vr1 Vr(1:N,j) Vr(N+1:2*N,j)];
end

%% orthogonalize the split sets using Grahm-Schmidt
Bc = my_gs(Vc1,10*eps);
Br = my_gs(Vr1,10*eps);

nc = size(Bc,2);
nr = size(Br,2);

Bc_ortho_diff = norm(eye(nc) - Bc'*Bc,inf);
Br_ortho_diff = norm(eye(nr) - Br'*Br,inf);

%  Just in case, we can compare with orthogonalization via qr
% [Bcqr,Rcqr] = qr(Vc1,0);
% [Brqr,Rrqr] = qr(Vr1,0);


%% Negate some columns of bases
%  to make sure that corresponding columns are comparable, we check the sign
%  of the largest entry (in absolute value) in each column of Bc and Br,
%  and replace the column by its negative if the sign is negative.
for j = 1:nc,
    [dummy,i] = max(abs(Bc(:,j)));
    if Bc(i,j) < 0,
        Bc(:,j) = - Bc(:,j);
    end
end

for j = 1:nr,
    [dummy,i] = max(abs(Br(:,j)));
    if Br(i,j) < 0,
        Br(:,j) = - Br(:,j);
    end
end

d = norm(Bc - Br,inf);


%% ------------------------------------------------------------
    function y = MultA(x)
        %  Computes y = A * x
        %  Roland W. Freund
        %  last change:  April 14, 2005
        
        y = A * x;
        
    end

    function y = MultAr(x)
        %  Computes y = Ar * x
        %  Roland W. Freund
        %  last change:  April 14, 2005
        
        y = Ar * x;
    end

    function y = MultAr2(x)
        %  Computes y = Ar * x
        %     implicitly using complex matrix-vector product
        xc = complex(x(1:N), x(N+1:end));
        tmp = A * xc;
        y = [real(tmp); imag(tmp)];
    end

end