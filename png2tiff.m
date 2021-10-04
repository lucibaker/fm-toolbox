function png2tiff(fn1,fn2)
% convert fn1.png to fn2.tiff 

A = imread([fn1 '.png']);
imwrite(A,['../productionfiles/' fn2 '.tiff']);