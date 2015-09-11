function s0 = getpoint(ax,ROI)


axis(ax);
h = line([0 0],2*pi*10.^ROI,'lineWidth',1,'Color','r','LineStyle','-'); % segment on i-axis
putaxis;

[x y] = ginput(1);
s0 = complex(x,y);
