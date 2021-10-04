function matR = rotate_mat(theta)
% This function creates the homogeneous transformation matrix for a 2D
% rotation of theta radians.

ctheta = cos(theta);
stheta = sin(theta);

matR = [ctheta -stheta 0;
        stheta ctheta 0;
        0 0 1];
    