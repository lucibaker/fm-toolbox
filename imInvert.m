function imInvert(filename, no_of_frames)
% inverts images in a .img file and writes a new .img file
% filename: eg. run_a, run_b.cut.bsub
% no_of_frames: number of frames to invert

sigma = 10;
img = imgRead([filename '.img']);
img2 = imgNew(['inv_' filename '.img'],img);
parfor i = 0:no_of_frames-1
%     filestring = [filename '_' num2str(i,'%04d') '_Cam_9994_Cine1.bmp'];
%     A = imread(filestring);
    A = imgGetFrame(img,i);
    A = uint8(round(imgaussfilt(255-A,sigma)));
%     imwrite(A,['inv_' filestring],'bmp');
    imgPutFrame(img2,A);
end

end