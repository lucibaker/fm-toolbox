function [N,C] = jpdf(var1, var2, nbins, varargin)
% Returns and plots the joint pdf of two variables var1 and var2.

% Inputs:
% var1: nx1 vector containing first variable (x-axis)
% var2: nx1 vector containing second variable (y-axis)
% nbins: number of bins to use in the JPDF
% varargin{1}: 4x1 vector [xlims ylims] specifying bounds of var1 and var2
% (optional)
% varargin{2}: 1 for pcolor, 0 for contourf (default) (optional)
% varargin{3}: colormap (default: hot) (optional)

% Outputs:
% N: [nbins x nbins] matrix containing the JPDF observation counts
% C: [nbins x 2] matrix containing the centers of the JPDF bins

if ~isempty(varargin)
    % data limits
    lims = varargin{1};
    idx = var1 >= lims(1) & var1 <= lims(2) & var2 >= lims(3) & var2 <= lims(4);
    var1 = var1(logical(idx)); var2 = var2(logical(idx));
    % pcolor or contourf
    if length(varargin) > 1
        pcolor_flag = varargin{2};
    else
        pcolor_flag = 0;
    end
    % colormap
    if length(varargin) > 2
        cmap_option = varargin{3};
    else
        cmap_option = 'hot';
    end
else
    pcolor_flag = 0;
    cmap_option = 'hot';
end


% C = {logspace(log10(min(var1)),log10(max(var1)),nbins)',...
%     linspace(min(var2),max(var2),nbins)'};    % histogram edge vector
[N,C] = hist3([var1, var2], [nbins,nbins]);
% hist3([var1, var2], [nbins,nbins]);
N = N'/sum(N(:));
ncont = 20; % number of contours
% vcont = logspace(log10(10e-7),log10(0.1*max(N(:))),ncont); % min and max contour levels (optional)
vcont = linspace(0,max(N(:)),ncont); % min and max contour levels (optional)

if pcolor_flag
    pcolor(C{1,1},C{1,2},N); shading flat  % ncont contour divisions
else
    contourf(C{1,1},C{1,2},N,[vcont],'LineColor','none')  % ncont contour divisions
end
c=colormap(cmap_option);
c2=zeros(length(c),3);
for i=1:length(c)
c2(i,:) = c(end-i+1,:);
end
colormap(c2);

axis('tight')

end