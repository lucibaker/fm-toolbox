function subplot_labels(fh,RelX,RelY)
% add subplot labels (a, b, c) to subplots in figure with handle fh at positions RelX and RelY
% RelX = .02; RelY = .9;
labelstr = {'({\it a})','({\it b})','({\it c})','({\it d})','({\it e})','({\it f})',...
    '({\it g})','({\it h})','({\it i})','({\it j})','({\it k})','({\it l})','({\it m})'};
sp = findobj(fh,'Type','Axes'); sp = flipud(sp);
for i = 1:length(sp)
    subplot(sp(i));
    DataX = interp1( [0 1], xlim(), RelX );
    DataY = interp1( [0 1], ylim(), RelY );
    text(DataX, DataY, labelstr{i},'interp','latex','fontsize',9)
end
