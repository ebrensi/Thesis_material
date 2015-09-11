close all
for j=[0.25 .5 1 1.5 2]
    figure;
    [hp htf] = rmp('ex1841s_3',5,j*i*pi*1e10);
    fname = sprintf('ex1841s_3_m5s0%d',j*100);
    
%     saveas(hp,sprintf('%s_pol',fname),'fig');
%     saveas(htf,sprintf('%s_tf',fname),'fig');
%     
%     exportfig(hp,sprintf('%s_pol.png',fname));
%     exportfig(htf,sprintf('%s_tf.png',fname));
end

