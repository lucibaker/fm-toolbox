function status = mask_circle_area(centers,Rp,h,status)
% mask PIV vectors inside detected hydrogel spheres

x_vec = h.x0:h.dx:(h.ix-h.dx);   % px
y_vec = h.y0:h.dy:(h.iy-h.dy);   % px
[x_vec_grid, y_vec_grid] = meshgrid(x_vec,y_vec);
parfor i = 1:h.nt
    insphere = zeros(length(y_vec),length(x_vec));
    for n = 1:size(centers{i},1)
        insphere = insphere + (sqrt((centers{i}(n,1)-x_vec_grid).^2 ...
            + (centers{i}(n,2)-y_vec_grid).^2) <= 1.1*Rp);
    end
    s = status(:,:,i);
    s(logical(insphere)) = 5;
    status(:,:,i) = s;
end