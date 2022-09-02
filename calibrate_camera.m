function rectify = calibrate_camera(I,W,order)
% undistortion and rectification of points in an image 
%
% rectify: function handle to map image coordinates to world coordinates 
% I: set of known calibration points in image coordinates (n x 2 vector) [px]
% W: set of known calibration points in world coordinates (n x 2 vector) [m]
% order: 1 for linear transformation (corrects for camera viewing angle but not lens distortion), 
%        2 for quadratic transformation (corrects for camera viewing angle and lens distortion)
%
% references: Fujita et al 1998 (Water Res.), Creutin et al 2003 (J. Hydrol.)


% find transformation coefficients

if order == 2
    A = [I.^2, I, ones(size(I,1),1), -I(:,1).^2.*W(:,1), -I(:,2).^2.*W(:,1), -I(:,1).*W(:,1), -I(:,2).*W(:,1), zeros(size(I,1),5);
        zeros(size(I,1),5), -I(:,1).^2.*W(:,2), -I(:,2).^2.*W(:,2), -I(:,1).*W(:,2), -I(:,2).*W(:,2), I.^2, I, ones(size(I,1),1)];
else
    A = [I, ones(size(I,1),1), -I(:,1).*W(:,1), -I(:,2).*W(:,1), zeros(size(I,1),3);
        zeros(size(I,1),3), -I(:,1).*W(:,2), -I(:,2).*W(:,2), I, ones(size(I,1),1)];
end
Z = [W(:,1); W(:,2)];
B = (A'*A)^-1*A'*Z;
    

% function to map image coords to world coords

if order == 2
    rectify = @(I) [(B(1)*I(:,1).^2 + B(2)*I(:,2).^2 + B(3)*I(:,1) + B(4)*I(:,2) + B(5))./ ...
        (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1), ...
        (B(10)*I(:,1).^2 + B(11)*I(:,2).^2 + B(12)*I(:,1) + B(13)*I(:,2) + B(14))./ ...
        (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1)];

else
    rectify = @(I) [(B(1)*I(:,1) + B(2)*I(:,2) + B(3))./ ...
        (B(4)*I(:,1) + B(5)*I(:,2) + 1), ...
        (B(6)*I(:,1) + B(7)*I(:,2) + B(8))./ ...
        (B(4)*I(:,1) + B(5)*I(:,2) + 1)];
end

end