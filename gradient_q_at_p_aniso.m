 function [dqdx,dqdy,valid_frac] = gradient_q_at_p_aniso(xp,yp,r1,r2,theta_p,xf,yf,q,dmin,sigma,theta_cut)
% interpolate fluid quantity q at anisotropic particle location xp, yp 
% using method in Schneiders et al. 2019 (approximates particles as ellipsoids)
%
% xp, yp: particle coordinates (scalars)
% r1, r2: imaged length of particle semimajor and semiminor axes, respectively
% theta_p: imaged particle tilt (in-plane orientation of major axis w.r.t. horizontal axis) [rad]
% xf, yf: fluid coordinates (vectors)
% q: fluid quantity (matrix)
% dmin: offset of interpolation band from particle
% sigma: width of interp. band
% theta_cut: cutoff angle above and below horz. axis of downwind portion of interp. band [rad]
% dqdx, dqdy: gradient of q at particle location in x and y directions (scalar)
% valid_frac: fraction of vectors in interp. band that were valid

if all(isreal([r1 r2 theta_p]))

    % pad q
    xf = xf(:); yf = yf(:);
    dxf = diff(xf(1:2)); dyf = diff(yf(1:2));
    q = padarray(q,[ceil((dmin+sigma)/dyf) ceil((dmin+sigma)/dxf)],nan,'both');
    xf = [((xf(1)-ceil((dmin+sigma)/dxf)*dxf):dxf:(xf(1)-dxf))'; xf; ((xf(end)+dxf):dxf:(xf(end)+ceil((dmin+sigma)/dxf)*dxf))'];
    yf = [((yf(1)-ceil((dmin+sigma)/dyf)*dyf):dyf:(yf(1)-dyf))'; yf; ((yf(end)+dyf):dyf:(yf(end)+ceil((dmin+sigma)/dyf)*dyf))'];
    [xf,yf] = meshgrid(xf,yf);

    % get elliptical boundaries interpolation band
    bound1 = ((xf-xp)*cos(theta_p) + (yf-yp)*sin(theta_p)).^2/(r1+dmin)^2 + ...
        ((xf-xp)*sin(theta_p) + (yf-yp)*cos(theta_p)).^2/(r2+dmin)^2 - 1;
    bound2 = ((xf-xp)*cos(theta_p) + (yf-yp)*sin(theta_p)).^2/(r1+dmin+sigma)^2 + ...
        ((xf-xp)*sin(theta_p) + (yf-yp)*cos(theta_p)).^2/(r2+dmin+sigma)^2 - 1;
    bound3 = abs(atan2(yf-yp, xf-xp));
    idx = logical(bound1 > 0 & bound2 < 0 & bound3 > theta_cut);

    q_n = q;
    q_n(~idx) = nan;
    
    % take gradients within interpolation band
    [dqdx_n,dqdy_n] = gradient(q_n,dxf);
    dqdx = nanmean(dqdx_n(:));
    dqdy = nanmean(dqdy_n(:));

    % valid fraction
    valid_frac = sum(~isnan(q_n(:)))/sum(idx(:));
    if valid_frac == 0
        dqdx = nan;
        dqdy = nan;
    end

%     % plot
%     r = 4;
%     plot_lims = [xp-r*(r1+dmin+sigma), xp+r*(r1+dmin+sigma), yp-r*(r1+dmin+sigma), yp+r*(r1+dmin+sigma)];
%     plot_idx_x = xf(1,:) > plot_lims(1) & xf(1,:) < plot_lims(2);
%     plot_idx_y = yf(:,1) > plot_lims(3) & yf(:,1) < plot_lims(4);
%     contourf(xf(plot_idx_y,plot_idx_x),yf(plot_idx_y,plot_idx_x),q(plot_idx_y,plot_idx_x),20,'edgecolor','none'); colormap jet; hold on
%     quiver(xf,yf,q,zeros(size(q)))
%     ellipse(r1,r2,theta_p,xp,yp,'k');
%     ellipse(r1+dmin,r2+dmin,theta_p,xp,yp,'k');
%     ellipse(r1+dmin+sigma,r2+dmin+sigma,theta_p,xp,yp,'k');
%     line([xp xp+r1+dmin+sigma],[yp yp+tan(theta_cut)*(r1+dmin+sigma)],'color','k')
%     line([xp xp+r1+dmin+sigma],[yp yp-tan(theta_cut)*(r1+dmin+sigma)],'color','k')
%     axis equal; axis(plot_lims); 
%     c = colorbar; c.Label.String = 'u [m/s]';
%     goodplot([5 4]);
%     hold off
%     keyboard
    
else
    % in case of imaginary r1, r2, theta_p
    valid_frac = nan;
    dqdx = nan;
    dqdy = nan;
end

end