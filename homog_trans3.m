function T = homog_trans3(a, b, c, xt)
% returns the 3D homogeneous transformation matrix for the following
% transformations (rotation, then translation):
% 1. rotate about x axis by [a] radians CCW
% 2. rotate about y axis by [b] radians CCW
% 3. rotate about z axis by [c] radians CCW
% 4. translate by 3x1 vector [xt] = [xt yt zt] 

T = [cos(a)*cos(b), cos(a)*sin(b)*sin(c) - sin(a)*cos(c), cos(a)*sin(b)*cos(c) + sin(a)*sin(c), xt(1);
     sin(a)*cos(b), sin(a)*sin(b)*sin(c) + cos(a)*cos(c), sin(a)*sin(b)*cos(c) - cos(a)*sin(c), xt(2);
     -sin(b), cos(b)*sin(c), cos(b)*cos(c), xt(3);
     0, 0, 0, 1];