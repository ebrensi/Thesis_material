function hp = plotfullerr(ferr,tfunc_conv_tol)
hp = contourf(frq,nvals,log10(ferr).',32); % use contourf or pcolor
set(gca,'YDir','normal');
shading flat;
caxis([log10(tfunc_conv_tol) 0]);
colormap('hot')
colorbar;
title(txt)
end
