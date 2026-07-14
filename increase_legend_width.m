function increase_legend_width(lgd, width_factor)
%INCREASE_LEGEND_WIDTH Increase legend width.

lgd.Units = 'normalized';
pos = lgd.Position;
pos(3) = width_factor * pos(3);
lgd.Position = pos;

end