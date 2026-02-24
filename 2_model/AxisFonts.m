function [  ] = AxisFonts( xTickSize , xLabelSize , yTickSize , yLabelSize , LegendSize , TitleSize )

%%% sets font of plot

% x axis
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', xTickSize )
set(xl, 'FontSize', xLabelSize);

% y axis
yl = get(gca,'YLabel');
yAY = get(gca,'YAxis');
set(yAY,'FontSize', yTickSize )
set(yl, 'FontSize', yLabelSize);

% legend
set(legend,'fontsize',LegendSize);

% title
ax = gca;
ax.TitleFontSizeMultiplier = TitleSize / 10 ; % 10 is default fontsize

end

