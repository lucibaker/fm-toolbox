function [q_at_p,valid_frac] = interp_q_at_p(xp,yp,xf,yf,q,r,IDW_flag)
% interpolate fluid quantity q at particle location xp, yp using either
% inverse-distance weighting or simple moving average
%
% xp, yp: particle coordinates (scalars)
% xf, yf: fluid coordinates (vectors)
% q: fluid quantity (matrix)
% r: neighborhood size (eg., r=2 uses the 2x2 square surrounding the
%   particle to interpolate)
% IDW_flag: 1 for inverse-distance weighting, 0 for simple moving average
% q_at_p: q interpolated at particle location (scalar)
% valid_frac: fraction of vectors in neighborhood that were valid

% pad q
q = padarray(q,[ceil(r/2) ceil(r/2)],nan,'both');
xf = xf(:); yf = yf(:);
dxf = diff(xf); dyf = diff(yf);
xf = [((xf(1)-ceil(r/2)*dxf):dxf:(xf(1)-dxf))'; xf; ((xf(end)+dxf):dxf:(xf(end)+ceil(r/2)*dxf))'];
yf = [((yf(1)-ceil(r/2)*dyf):dyf:(yf(1)-dyf))'; yf; ((yf(end)+dyf):dyf:(yf(end)+ceil(r/2)*dyf))'];

% get neighborhood
[Bx,I] = sort(abs(xf-xp),'ascend');
x_idx = I(1:r);
[By,I] = sort(abs(yf-yp),'ascend');
y_idx = I(1:r);

q_n = q(y_idx,x_idx);

% keyboard
% simple moving average
if ~IDW_flag
    q_at_p = nanmean(q_n(:));

% inverse-distance weighting
else
    [dx,dy] = meshgrid(Bx(1:r),By(1:r));
    inv_dists = 1./sqrt(dx.^2 + dy.^2);
    total_inv_dist = sum(sum(inv_dists(~isnan(q_n))));
    q_at_p = nansum(q_n(:).*inv_dists(:))/total_inv_dist;    
end

% valid fraction
valid_frac = sum(~isnan(q_n(:)))/r^2;
if valid_frac == 0
    q_at_p = nan;
end

end