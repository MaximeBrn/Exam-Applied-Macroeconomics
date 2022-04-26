function f_plot_with_CI(x_axis,res_mean,res_CI_lower,res_CI_upper,str_title,interpreter,zero_line)
    t_conf = [x_axis fliplr(x_axis)];
    y_conf=[res_CI_lower' flipud(res_CI_upper)'];
    zeroline = ones(length(x_axis),1)*0;


    if zero_line==1
        plot(x_axis,res_mean,'LineWidth',1.5,'color','black'); hold on
        plot(x_axis,zeroline,'LineWidth',1); hold on
    else 
        plot(x_axis,res_mean,'LineWidth',1.5,'color','black'); hold on
    end
    
    fill(t_conf, y_conf, 1,'facecolor', 'black', 'edgecolor', 'black', 'facealpha', 0.3,'edgealpha',0.8);
    xlabel('quarters','FontSize',8);
    if interpreter=="Latex"
        title(str_title,'interpreter','latex','FontSize',10); 
    else
        title(str_title,'FontSize',10);
    end


end
