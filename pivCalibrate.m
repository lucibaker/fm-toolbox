function pxtom = pivCalibrate()
% Calculate pxtom from calibration image
% takes the filename of the calibration bmp image and returns the scale
% factor in m/pixel.

% Inputs
numrows = 23; %23; % 12; %   % number of rows of dots
numcols = 18;     % number of columns of dots
rrange = [6 9];  %[8 12]; % %radius range of dots [px]
% threshold = 100; %80; %  background brightness of image (0 to 255)
sens = 0.96;    % imfindcircles sensitivity level
fn = '../calibration.bmp';     % filename of calibration image

% Process calibration image
A = imread(fn);
% A(A<threshold)=0;
% A(A>threshold)=255;
warning('off','images:initSize:adjustingMag');
[cc,rr] = imfindcircles(A,rrange,'Sensitivity',sens);
figure; imshow(A); title('Calibration')
viscircles(cc,rr,'Edgecolor','b');

% Compute mean distance between dots
if numrows*numcols ~= length(rr)
    error('Number of dots detected does not equal numrows x numcols.')
end
rows = sortrows(cc,2); rows = rows(floor(numrows/2)*numcols+1:(floor(numrows/2)+1)*numcols,1); 
rows = diff(sortrows(rows));
cols = sortrows(cc,1); cols = cols((floor(numcols/2)*numrows+1):(floor(numcols/2)+1)*numrows,2); 
cols = diff(sortrows(cols));
pxtom = vertcat(rows,cols);

% scale factor [m/px]
pxtom = (2.5e-3)/mean(pxtom);

end