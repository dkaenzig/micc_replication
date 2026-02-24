function [] = plotCI(x,y1,y2,color,alpha)

%%% helper that plots confidence bands

fill([x', fliplr(x')], [y1, fliplr(y2)], color, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end