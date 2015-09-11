function v = getview(ax)

if nargin == 0
	ax = gca;
end
v.pba = get(ax, 'PlotBoxAspectRatio');
v.dar = get(ax, 'DataAspectRatio');
v.cva = get(ax, 'CameraViewAngle');
v.cuv = get(ax, 'CameraUpVector');
v.ct = get(ax, 'CameraTarget');
v.cp = get(ax, 'CameraPosition');
