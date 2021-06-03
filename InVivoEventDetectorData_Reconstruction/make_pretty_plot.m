function[h_fig_out] = make_pretty_plot(t,ymat,h_fig_in,colorvec,ylimits,xlimits)

meany = mean(ymat,2);
errory = std(ymat,0,2);
 
h_fig_out = figure(h_fig_in);
hold on        
h_ax = errorbar(t,meany,errory,'.--','MarkerSize',30);
set(h_ax,'MarkerFaceColor',colorvec,'MarkerEdgeColor',colorvec)
set(h_ax,'Color',[ .4 .4 .4]);
    
set(gca,'FontSize',40);
ylim(ylimits);
xlim(xlimits)

hold off