function matT = translate_mat(dx,dy)
% This function creates the homogeneous transformation matrix for a 2D
% translation by dx in the x direction and dy in the y direction.

matT = [1 0 dx;
        0 1 dy;
        0 0 1];