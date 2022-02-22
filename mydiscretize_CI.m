function [xx,yy]=mydiscretize_CI(xdat,ydat,func,varargin)
% bins input data and outputs expected values and 95% confidence intervals
% for the input function, using bootstrapping

%can input number of bins (N), which will be based on bins with constant
%sample sizes
%mydiscretize_CI(xdat,ydat,func,N)

%can also input bin edges/limits directly:
%mydiscretize_CI(xdat,ydat,func,[],edges)

%plot output xx,yy with function errorbar_CI

minbinsize=3; %3

if ~isempty(varargin)
    N=varargin{1};
else
    N=10;
end
lims=[];
if length(varargin)>1
    lims=varargin{2};
end
if isempty(lims)

sorted=sort(xdat);
sorted=sorted(~isnan(sorted));
ind=round(linspace(1,length(sorted),N));
lims=sorted(ind);
%lims=linspace(sorted(1),sorted(end),N);
end

if length(varargin)>2
    minbinsize=varargin{3};
end
[Y,~]=discretize(xdat,lims);

bins=unique(Y);
bins=bins(~isnan(bins));

xx=nan(max(bins),3);
yy=xx;

for i=bins(:)'

    ind=i==Y;
    
    if sum(ind)<minbinsize
        continue
    end
    try
    [CI avg]=bootstrap_CI(ydat(ind),func);
    yy(i,:)=[avg CI];
    catch
        keyboard
    end

    
    
    [CI avg]=bootstrap_CI(xdat(ind),@mean);
    xx(i,:)=[avg CI];

  
end

end
