function putaxis(ax)

if ~nargin
	ax = axis();
end

xmin = ax(1);
xmax = ax(2);
ymin = ax(3);
ymax = ax(4);

h = line([xmin xmax],[0 0],'lineWidth',1,'Color','g','LineStyle','--');
uistack(h,'bottom');
h = line([0 0],[ymin ymax],'lineWidth',1,'Color','g','LineStyle','--');
uistack(h,'bottom');

