function img_to_bmp(run_id)
% convert .img file to series of .bmp images
addpath('Z:\luci\MATLAB\pivio-matlab');

% run_id = 'run-d';
img = imgRead([run_id '.bsub.cut.fill.img']);
parfor i = 1:img.it
    A = uint8(imgGetFrame(img,i-1));
    imwrite(A,['images_masked/' run_id sprintf('-%04d',i) '.bmp']);
end