function [a_cond, b_centers, c_centers, N] = condition_vars2(a,b,c,nbins,scaleflag,varargin)
% condition a on b and c 
% divide b and c into nbins bins (2x1 vector)
% scaleflag: 0=lin, 1=log spaced bins (2x1 vector)
% blims: (optional) specify the min and max b, min and max c lims (4x1
% vector)
% a_cond: nbins(1) x nbins(2) matrix
% N is the number of data points in each bin

if ~isempty(varargin)
    bclims = varargin{1};
else
    bclims = [min(b),max(b),min(c),max(c)];
end

if scaleflag(1)
    b_bins = logspace(log10(bclims(1)),log10(bclims(2)),nbins(1)+1);
else
    b_bins = linspace(bclims(1),bclims(2),nbins(1)+1);
end
if scaleflag(2)
    c_bins = logspace(log10(bclims(3)),log10(bclims(4)),nbins(2)+1);
else
    c_bins = linspace(bclims(3),bclims(4),nbins(2)+1);
end
a_cond = zeros(nbins(1),nbins(2));
N = zeros(nbins(1),nbins(2));
for n = 1:nbins(1)
    for m = 1:nbins(2)
        idx = b>b_bins(n) & b<b_bins(n+1) & c>c_bins(m) & c<c_bins(m+1);
        a_cond(n,m) = nanmean(a(idx));
        N(n,m) = sum(idx);
    end
end
b_centers = b_bins(2:end)-diff(b_bins)/2;
c_centers = c_bins(2:end)-diff(c_bins)/2;
end