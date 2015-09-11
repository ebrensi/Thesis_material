function [err fullerr] = tfunc_err(URM_FR,ROM_FR)
		Ls = size(URM_FR,3);
		fullerr = zeros(1,Ls);
		for i = 1:Ls
			fullerr(i) = norm(URM_FR(:,:,i) - ROM_FR(:,:,i)) / norm(URM_FR(:,:,i));
		end
		err = norm(fullerr);
end

	
