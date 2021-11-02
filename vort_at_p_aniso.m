 function [vort_at_p,valid_frac] = vort_at_p_aniso(xp,yp,r1,r2,theta_p,xf,yf,u,v,dmin)
% calculate vorticity of fluid velocity at anisotropic particle located at xp, yp using Stokes' theorem 
%
% xp, yp: particle coordinates (scalars)
% r1, r2: imaged length of particle semimajor and semiminor axes, respectively (scalars)
% theta_p: imaged particle tilt (in-plane orientation of major axis w.r.t. horizontal axis) [rad] (scalars)
% xf, yf: fluid coordinates (vectors)
% u: fluid u velocity (matrix)
% v: fluid v velocity (matrix)
% dmin: offset of integration curve from particle
% vort_at_p: fluid vorticity at particle location (scalar)
% valid_frac: fraction of vectors in interp. band that were valid

if all(isreal([r1 r2 theta_p]))
    % get image bounds
    img_bounds = [xf(1) xf(end) yf(1) yf(end)];
    
    % pad u, v
    xf = xf(:); yf = yf(:);
    dxf = diff(xf(1:2)); dyf = diff(yf(1:2));
    u = padarray(u,[ceil(dmin/dyf) ceil(dmin/dxf)],nan,'both');
    v = padarray(v,[ceil(dmin/dyf) ceil(dmin/dxf)],nan,'both');
    xf = [((xf(1)-ceil(dmin/dxf)*dxf):dxf:(xf(1)-dxf))'; xf; ((xf(end)+dxf):dxf:(xf(end)+ceil(dmin/dxf)*dxf))'];
    yf = [((yf(1)-ceil(dmin/dyf)*dyf):dyf:(yf(1)-dyf))'; yf; ((yf(end)+dyf):dyf:(yf(end)+ceil(dmin/dyf)*dyf))'];
    [xf,yf] = meshgrid(xf,yf);

%     % get elliptical integration curve 
%     t = linspace(0,2*pi,50);
%     int_curve_x = xp + (r1+dmin)*cos(theta_p)*cos(t) - (r2+dmin)*sin(theta_p)*sin(t);
%     int_curve_y = yp + (r1+dmin)*sin(theta_p)*cos(t) + (r2+dmin)*cos(theta_p)*sin(t);
%     
% %     ((xf-xp)*cos(theta_p) + (yf-yp)*sin(theta_p)).^2/(r1+dmin)^2 + ...
% %         ((xf-xp)*sin(theta_p) + (yf-yp)*cos(theta_p)).^2/(r2+dmin)^2 - 1;
%     
%     u_curve = interp2(xf,yf,u,int_curve_x,int_curve_y);
%     v_curve = interp2(xf,yf,v,int_curve_x,int_curve_y);
%     
%     dl_x = ?;
%     dl_y = ?;
    
    % get rectangular integration curve
    x_bound1 = xp - abs((r1+dmin)*cos(theta_p));
    x_bound2 = xp + abs((r1+dmin)*cos(theta_p));
    y_bound1 = yp - abs((r1+dmin)*sin(theta_p));
    y_bound2 = yp + abs((r1+dmin)*sin(theta_p));
    
    [~,x1_idx] = min(abs(xf(1,:) - x_bound1)); 
    [~,x2_idx] = min(abs(xf(1,:) - x_bound2));
    [~,y1_idx] = min(abs(yf(:,1) - y_bound1)); 
    [~,y2_idx] = min(abs(yf(:,1) - y_bound2));
    
    u_curve = [u(y1_idx,x1_idx:x2_idx)'; -u(y2_idx,x1_idx:x2_idx)'];
    v_curve = [-v(y1_idx:y2_idx,x1_idx); v(y1_idx:y2_idx,x2_idx)];
    u_curve = [zeros(size(v_curve)); u_curve]; v_curve = [v_curve; zeros(size(u_curve))];
    xf_curve = [xf(y1_idx:y2_idx,x1_idx); xf(y1_idx:y2_idx,x2_idx); xf(y1_idx,x1_idx:x2_idx)'; xf(y2_idx,x1_idx:x2_idx)'];
    yf_curve = [yf(y1_idx:y2_idx,x1_idx); yf(y1_idx:y2_idx,x2_idx); yf(y1_idx,x1_idx:x2_idx)'; yf(y2_idx,x1_idx:x2_idx)'];
    
    % check for circulation integrals that go out of bounds
    if any(xf_curve < img_bounds(1)+dxf | xf_curve > img_bounds(2)-dxf | yf_curve < img_bounds(3)+dyf | yf_curve > img_bounds(4)-dyf)
        valid_frac = nan;
        vort_at_p = nan;   
    else
    
        % interpolate NVVs
        valid_frac = numel(xf_curve);
        for i = 1:length(u_curve)
            if isnan(u_curve(i))
                valid_frac = valid_frac - 1;
                interp_idx = [find(xf(1,:)==xf_curve(i)), find(yf(:,1)==yf_curve(i))];
                u_curve(i) = nanmean(u(interp_idx(2)-1:interp_idx(2)+1,interp_idx(1)-1:interp_idx(1)+1),'all');
                v_curve(i) = nanmean(v(interp_idx(2)-1:interp_idx(2)+1,interp_idx(1)-1:interp_idx(1)+1),'all');
                if i < y2_idx-y1_idx+1
                    v_curve(i) = -v_curve(i);
                else
                    u_curve(i) = -u_curve(i);
                end
            end
        end
        valid_frac = valid_frac/numel(xf_curve);

        % get circulation (line integral of fluid vel around curve)
        circ = dyf*sum(v_curve) + dxf*sum(u_curve);

    %     dyf*(-sum(v(y1_idx:y2_idx,x1_idx)) + sum(v(y1_idx:y2_idx,x2_idx))) + ...
    %         dxf*(sum(u(y1_idx,x1_idx:x2_idx)) - sum(u(y2_idx,x1_idx:x2_idx)));

        vort_at_p = circ/((x_bound2-x_bound1)*(y_bound2-y_bound1)); % divide circ by area enclosed by integration curve
    
    end
    
else
    % in case of imaginary r1, r2, theta_p
    valid_frac = nan;
    vort_at_p = nan;
end

end