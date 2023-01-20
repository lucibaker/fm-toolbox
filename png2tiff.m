% convert png's to tiff's in current folder
D = dir;
D = D(3:end);

for i = 1:length(D)
    fn = D(i).name(1:end-4);
    png2tiff_convert(fn,fn);
end

function png2tiff_convert(fn1,fn2)
% convert fn1.png to fn2.tiff 
A = imread([fn1 '.png']);
imwrite(A,['../figures-tiff/' fn2 '.tiff']);
end