function multi_ex(k,todo)
% perform some analysis on several example data sets

% warning off
examples = {'1a','1b','2','308_1','308_2','308_3',...
           '1346_1','1346_2','1346_3','1346_4',...
           '13875_1','13875_2','13875_3',...
            '6','19','30','240a','240b','448'};

example =  examples(todo);
n_e = length(example);

% err_hess = zeros(k,n_e);
% err_proj = zeros(k,n_e);
for e = 1:n_e
%     example{e}
    data = sprintf('ex%s',char(example(e)));

%% which examples have b==c ? 
% [A E c b] = inputdata(data);
% if all(b == c)
%     fprintf('%s\n',data);
% end

% 1a, 
% 308_1, 308_3, 
% 1346_1, 1346_4, 
% 13875_1, 6, 
% 19, 30, 240a, 240b, 448 (all lossless positive real examples)
  
%% k vs s0
%     [min_K s0_opt] = k_vs_s0_1D(data,100);
%     saveas(gcf,sprintf('%s_k_vs_s0_1D',data));
% k_vs_s0_2D(data,k); 

%% tfunc err vs k
%     err = err_vs_k(data,k);   
%     err_hess(:,e) = err(:,1);
%     err_proj(:,e) = err(:,2);
    
%% compute transfer function 
%     [C G c b] = inputdata(data);
%     tfunc = abs(tfunc_urm(data,C,G,c,b));
%     frq = logspace(8,10,200); % 200 points from 10^8 to 10^10, distributed logly
%     plot_tfunc(frq,tfunc);
%     title(data);
%     saveas(gcf,sprintf('%s_tfunc_urm',data));

%%  plot err vs s0 (real s0)
% err_vs_s0_1D(data,100,100);


%%  plot err vs s0 (complex s0)
%     err_vs_s0_2D(data,k,50,'hess');
%     err_vs_s0_2D(data,k,50,'proj'); % tfunc error

%% plot cond(Ck-sGk) vs s vs k
%     title(sprintf('%s:cond(C_k-sG_k)',data))
%     saveas(gcf,sprintf('condH%s.fig',data));
%     drawnow

end

% for use with err vs k
% semilogy(1:k,err_hess);
% title(sprintf('hess err vs k'));
% legend(example);
% figure;
% semilogy(1:k,err_proj);
% title(sprintf('PRIMA err vs k'));
% legend(example);
% --


% semilogy(1:k,cCs);
% title(sprintf('cond(C_k)'));
% legend(example);
% 
% figure;
% semilogy(1:k,cGs);
% title(sprintf('cond(G_k)'));
% legend(example)

end
