function [dqdx,dqdy,valid_frac] = gradient_q_at_p(xp,yp,xf,yf,q,r)
% calculate gradient of fluid quantity q at particle location xp, yp using
% central differences in a rxr neighborhood surrounding particle. The mean
% is then taken over the neighborhood. Fluid gridpoints must be equally
% spaced with the same spacing in x and y directions.
%
% xp, yp: particle coordinates (scalars)
% xf, yf: fluid coordinates (vectors)
% q: fluid quantity (matrix)
% r: neighborhood size (eg., r=2 uses the 2x2 square surrounding the
%   particle to interpolate)
% dqdx, dqdy: gradient of q at particle location in x and y directions 
%   (scalar)
% valid_frac: fraction of vectors in neighborhood that were valid

% get neighborhood
[~,I] = sort(abs(xf-xp),'ascend');
x_idx = I(1:r);
[~,I] = sort(abs(yf-yp),'ascend');
y_idx = I(1:r);

q_n = q(y_idx,x_idx);

% take gradient
dx = xf(2)-xf(1);
[dqdx_n,dqdy_n] = gradient(q_n,dx);

% simple moving average
dqdx = nanmean(dqdx_n(:));
dqdy = nanmean(dqdy_n(:));

% valid fraction
valid_frac = sum(~isnan(q_n(:)))/r^2;
if valid_frac == 0
    dqdx = nan;
    dqdy = nan;
end

end