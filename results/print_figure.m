function print_figure(handler_figure,name_file,x_size,y_size)
% set(get(handler_figure,'CurrentAxes'), 'Position', get(get(handler_figure,'CurrentAxes'), 'OuterPosition') - ...
%     get(get(handler_figure,'CurrentAxes'), 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
set(get(handler_figure,'CurrentAxes'),'LooseInset',get(get(handler_figure,'CurrentAxes'),'TightInset'))

set(handler_figure, 'PaperUnits', 'centimeters');
set(handler_figure, 'PaperSize', [x_size y_size]);
set(handler_figure, 'PaperPositionMode', 'manual');
set(handler_figure, 'PaperPosition', [0 0 x_size y_size]);
set(handler_figure, 'renderer', 'OpenGL');
print(handler_figure, '-dpdf', sprintf('%s.pdf', name_file));
print(handler_figure, '-dpng', sprintf('%s.png', name_file));
print(handler_figure, '-depsc2',sprintf('%s.eps', name_file));
end