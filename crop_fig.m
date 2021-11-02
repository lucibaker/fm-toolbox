function crop_fig(fn)
% crop figure and save as png

f = openfig([fn '.fig']);
pause(1/10)

goodplot([3.5 2])
ax = gca;

% set(ax, 'Position', [0.02,0.1,1,0.82])
% set(gcf,'OuterPosition',[50 50 680 700]);

% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

% [left bottom ax_width ax_height]
% ax.Position=[-.1 -.1 1.2 1.35];
% ax.Position=[-.05 -.05 1.15 1.2];
ax.Position=[-.05 -.05 1.15 1.2];

% print([fn '_test'],'-dpng','-r150')
print(fn,'-dpng','-r600')

% savefig(fn);
close(f)

end