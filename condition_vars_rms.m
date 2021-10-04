function [a_cond, b_centers, N] = condition_vars_rms(a,b,nbins,scaleflag,varargin)
% condition a on b, dividing b into nbins bins into lin or log bins as
% specified by scaleflag: 0=lin, 1=log, from blims(1) to blims(2)
% blims is an optional 2x1 vector specifying the min and max b
% N is the number of data points in each bin

if ~isempty(varargin)
    blims = varargin{1};
else
    blims = [min(b),max(b)];
end

if scaleflag
    b_bins = logspace(log10(blims(1)),log10(blims(2)),nbins+1);
else
    b_bins = linspace(blims(1),blims(2),nbins+1);
end
a_cond = zeros(nbins,1);
N = zeros(nbins,1);
for n = 1:nbins
    a_cond(n,:) = rms(a(b>b_bins(n) & b<b_bins(n+1)),'omitnan');
    N(n) = length(find(b>b_bins(n) & b<b_bins(n+1)));
end
b_centers = b_bins(2:end)-diff(b_bins)/2;
end