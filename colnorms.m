function nrms = colnorms(A,type)

n = size(A,2);
nrms = zeros(1,n);
if exist('type','var')
	for i = 1:n
		nrms(i) = norm(A(:,i),type);
	end
else
	for i = 1:n
		nrms(i) = norm(A(:,i));
	end
end
