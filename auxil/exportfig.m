function exportfig(h,figname,ext)
% exportfig(h,figname,ext)
%
% save a figure image to a file, without the white border

if isinf(h)
    h=sort(findall(0,'type','figure'));
end

for j = 1:length(h)
    figure(h(j));
    ax = get(h(j),'CurrentAxes');
    TI = get(ax,'TightInset');
    OP = get(ax,'OuterPosition');
    Pos = OP + [ TI(1:2), -TI(1:2)-TI(3:4) ];
    set(ax,'Position',Pos);
    if length(h) == 1
        fprintf('\\putfig{.50}{%s.%s}\n',figname,ext);
        saveas(h,figname,ext);
    else
        fname = sprintf('%s_%d',figname,j);
        fprintf('\\putfig{.50}{%s.%s}\n',fname,ext);
        saveas(h(j),fname,ext);
    end
end
