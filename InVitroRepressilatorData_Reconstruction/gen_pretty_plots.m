function[fig_hand1,fig_hand2] = gen_pretty_plots(data_obj,tfmodel,t,color,ind)

dt = data_obj.Ts;
U = data_obj.u;
Y = data_obj.y;

trmin = 25;
trmax = numel(t);

ms = 25;
noutputs = size(Y,2);
ninputs = size(U,2);

[Ypred,fit,x0] = compare(data_obj,tfmodel);
Ypred = Ypred.Y(:,1:noutputs);



fig_hand1 = figure();
hold on

for j = 1:noutputs
    plot(t(trmin:trmax),Y(trmin:trmax,j),'.','MarkerSize',ms,'Color',color{j})
    plot(t(trmin:trmax),Ypred(trmin:trmax,j),'-','LineWidth',1,'Color',color{j})
    ylim([-.35 .35])
    xlim([t(trmin) t(trmax)])
end
hold off

set(gca,'FontSize',40)
%xlabel('Time (minutes)')
%ylabel('Normalized Output (a.f.u.)')

% fig_hand2 = figure()
% hold on
% for j = 1:2
%     if j<ind
%        colorvec = color{j};
%        plot(t(trmin:trmax),U(trmin:trmax,j),'.','MarkerSize',ms,'Color',colorvec);
%     end
%     if j==ind && j~=3
%        colorvec = color{j+1};
%        plot(t(trmin:trmax),U(trmin:trmax,j),'.','MarkerSize',ms,'Color',colorvec);
%     end
%     if j>ind
%        colorvec = color{j+1};
%        plot(t(trmin:trmax),U(trmin:trmax,j),'.','MarkerSize',ms,'Color',colorvec);
%     end
% end
% 
% colorvec = [0 .5 0];
% plot(t(trmin:trmax),U(trmin:trmax,3)/10,'-','LineWidth',4,'Color',[0 0 0]);
% 
% 
% set(gca,'FontSize',40)
% %ylim([-1 1])
% 
% xlim([t(trmin) t(trmax)])
% 
% %xlabel('Time (minutes)')
% %ylabel('Normalized Input (a.f.u.)')
% hold off
    
