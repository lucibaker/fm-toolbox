function pxtom = pivCalibrate(d_m, numrows, numcols, rrange, sens, fn)
% Calculate pxtom from calibration image. Detects dots in the calibration
% image and returns the scale factor in m/pixel.
%
% % Inputs
% d_m = 2.5e-3;                 % spacing between calibration dots [m]
% numrows = 23;                 % number of rows of dots
% numcols = 18;                 % number of columns of dots
% rrange = [6 9];               % radius range of dots [px]
% sens = 0.96;                  % imfindcircles sensitivity level
% fn = 'calibration.bmp';       % filename and path of calibration image

% Find dots in calibration image
A = imread(fn);
[cc,rr] = imfindcircles(A,rrange,'Sensitivity',sens);

% View image and detected dots
warning('off','images:initSize:adjustingMag');
figure; imshow(A); title('Calibration')
viscircles(cc,rr,'Edgecolor','b');

% check that the right number of dots are detected
if numrows*numcols ~= length(rr)
    error('Number of dots detected does not equal numrows x numcols.')
end

% sort dots into rows and columns
rows = sortrows(cc,2); rows = rows(floor(numrows/2)*numcols+1:(floor(numrows/2)+1)*numcols,1); 
cols = sortrows(cc,1); cols = cols((floor(numcols/2)*numrows+1):(floor(numcols/2)+1)*numrows,2); 

% Compute mean distance between dots in pixels
d_rows = diff(sortrows(rows));
d_cols = diff(sortrows(cols));
d_px = vertcat(d_rows,d_cols);  % spacing between dots [px]

% scale factor [m/px]
pxtom = d_m/mean(d_px);

end