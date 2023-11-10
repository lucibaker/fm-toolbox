function l = compare_plots(u, y, c, mk, ms, ls, w_u, w_y, c_w_u, c_w_y)
% compare N different plots on the current axes.
%
% inputs:
% u: variable on horizontal axis (Nx1 cell)
% y: variable on vertical axis (Nx1 cell)
% c: color (Nx1 cell)
% mk: marker (Nx1 cell)
% ms: marker size (Nx1 array)
% ls: line style (Nx1 cell)
% w_u (optional): error on u for error bars (Nx1 cell) 
% w_y (optional): error on y for error bars (Nx1 cell) 
% c_w_u (optional): error bar color for u (Nx1 cell)
% c_w_y (optional): error bar color for y (Nx1 cell)
%   Note: all 4 error arguments must be specified, or none. If ebars should   
%   only be plotted on one variable, set the cell entries for other to 'NaN'.
% 
% outputs:
% l: handle to each line (Nx1 array)

% number of lines to plot
N = length(u);
l = zeros(N,1);

% check inputs
switch nargin
    case 6
        w_u = num2cell(nan(N,1));
        w_y = num2cell(nan(N,1));
    case 10
    otherwise
        error('Incorrect number of input arguments.')
end

% plot error bars
for i = 1:N
    if any(~isnan(w_u{i}))
        patch = fill([u{i}(:)+w_u{i}(:); flipud(u{i}(:)-w_u{i}(:))], [y{i}(:); flipud(y{i}(:))], c_w_u{i});
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.15);
%         errorbar(u{i},y{i},w_u{i},'horizontal','.','color',c_w_u{i},'CapSize',2,'linewidth',0.5); 
        hold on
    end
    if any(~isnan(w_y{i}))
        patch = fill([u{i}(:); flipud(u{i}(:))], [y{i}(:)+w_y{i}(:); flipud(y{i}(:)-w_y{i}(:))], c_w_y{i});
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.15);
%         errorbar(u{i},y{i},w_y{i},'vertical','.','color',c_w_y{i},'CapSize',2,'linewidth',0.5); 
        hold on
    end
end

% plot lines
for i = 1:N
    l(i) = plot(u{i},y{i},'color',c{i},'marker',mk{i},'markersize',ms(i),'linestyle',ls{i},'markerfacecolor',c{i},'linewidth',1.5); hold on
end

end
