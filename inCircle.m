function B = inCircle(centers,radius,image_width,image_height)
% finds which pixels in an image are inside a circle 
% inputs:
% circle centers [px], circle radius [px], image matrix height and width
% [px]
% output:
% matrix B of same size as A with a value for each pixel: 
% 1 if inside, 0 if outside the circle

[x_grid, y_grid] = meshgrid(1:image_width,1:image_height);
B = zeros(size(x_grid));
for i = 1:size(centers,1)
    B(sqrt((centers(i,1)-x_grid).^2 + (centers(i,2)-y_grid).^2) <= radius) = 1;
end
B = logical(B);
end