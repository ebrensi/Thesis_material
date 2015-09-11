function [H H_string H_latex]  = spyH(H,nrm_H)

n = length(H);

if ~exist('nrm_H','var')
	nrm_H = norm(H);
end

relmag = abs(H) / nrm_H;
small = relmag < sqrt(eps);
nearzero = relmag <= 1e-12;
H(small) = 2;
H(nearzero) = 0;
H(~small) = 1;

H_string = [];
H_latex = sprintf('\\begin{tabular}{r|%s}\n','c'*ones(1,n));
for i = 1:n
	H_latex = sprintf('%s&%d',H_latex,i);
end
H_latex = sprintf('%s\\\\\n',H_latex);
for i = 1:size(H,1)
	H_string = sprintf('%s%2d  ',H_string, i);
	H_latex = sprintf('%s%2d ',H_latex, i);
	for j = 1:size(H,2)
		val = H(i,j);
		if val == 0
			hval = ' ';
		elseif val == 1
			hval = '*';
		elseif val == 2
			hval = '.';
		end
		H_string = sprintf('%s  %s',H_string, hval);
		H_latex = sprintf('%s&%s', H_latex, hval);
	end
	H_string = sprintf('%s\n', H_string);
	H_latex = sprintf('%s\\\\\n', H_latex);
end

end
