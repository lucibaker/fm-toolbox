function goodplot(papersize)
% function which produces a nice-looking plot
% and sets up the page for nice printing
f = gcf;

fontsize = 9; 
margin = 0.5;
intrp = 'latex';  % 'tex'; % 
% fnt = 'Times';

% axes
h = findobj(f,'Type','Axes');
for i = 1:length(h)
    set(h(i),'FontSize',fontsize);
    set(get(h(i),'xlabel'),'FontSize', fontsize, 'Interpreter', intrp); %, 'FontName', fnt);
    set(get(h(i),'ylabel'),'FontSize', fontsize, 'Interpreter', intrp); %, 'FontName', fnt);
    set(get(h(i),'zlabel'),'FontSize', fontsize, 'Interpreter', intrp); %, 'FontName', fnt);
    set(get(h(i),'title'),'FontSize', fontsize, 'Interpreter', intrp); %, 'FontName', fnt);
    set(h(i),'ticklength',2*[0.0100 0.0250]);
    set(h(i),'TickLabelInterpreter',intrp); %, 'FontName', fnt)
    set(h(i),'Box','on');
    %set(h(i),'ActivePositionProperty','outerposition');
end

% legend
l = findobj(f,'Type','Legend');
if ~isempty(l)
    set(l,'FontSize',fontsize, 'Interpreter',intrp); %, 'FontName', fnt);
%     legendmarkeradjust(20);
end

% annotations
a = findall(f,'Tag','scribeOverlay');
if ~isempty(a)
    ac = get(a,'Children');
    for i = 1:length(ac)
        try
            set(ac{i},'FontSize',fontsize, 'Interpreter',intrp); %, 'FontName', fnt);
        catch
        end
    end
end 

% colorbars
c = findall(f,'Type','ColorBar');
if ~isempty(c)
    for i = 1:length(c)
        set(c(i),'FontSize',fontsize,'TickLabelInterpreter',intrp); %, 'FontName', fnt);
        set(c(i).Label,'Interpreter',intrp,'FontSize',fontsize);
    end
end

% figure
set(f,'color','w');
set(f,'PaperUnits','inches');
set(f,'PaperSize', papersize);
set(f,'PaperPosition',...
    [margin margin papersize(1)-2*margin papersize(2)-2*margin]);