function f_plot_response(x_axis,res_mean,res_CI_lower,res_CI_upper,str_title)
    t_conf = [x_axis fliplr(x_axis)];
    y_conf=[res_CI_lower' flipud(res_CI_upper)'];
    plot(x_axis,res_mean,'LineWidth',2); hold on
    fill(t_conf, y_conf, 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4);
    xlabel('quarters','FontSize',8);
    title(str_title,'FontSize',10); 
end
