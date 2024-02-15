function figure_out = plot_graphics(x, y, inp_title, inp_xlabel, inp_ylabel)
    figure_out = figure('DefaultAxesFontSize',14);
    plot(x,y,'r-o','MarkerSize',5,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r');
    grid on;
    xlabel(inp_xlabel, 'FontSize', 24);
    ylabel(inp_ylabel, 'FontSize', 24);
    title(inp_title, 'FontSize', 24);
    saveas(figure_out, ['graphics/' inp_title '.png'])
    close
end
