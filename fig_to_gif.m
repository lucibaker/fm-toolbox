function fig_to_gif(gif_fn,delay,res)
% write current figure as a frame of a gif
% gif_fn: .gif filename (including extension)
% delay: interval between frames [s] (1/framerate)
% res: resolution flag string, e.g. '-r450' for 450 dpi
 
if nargin < 3
    res = '-r300';
end
cdata = print('-RGBImage',res);
frame.cdata = cdata; frame.colormap = [];
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if ~exist(gif_fn,'file')
    imwrite(imind,cm,gif_fn,'gif','Loopcount',inf,'DelayTime',delay);
else
    imwrite(imind,cm,gif_fn,'gif','WriteMode','append','DelayTime',delay);
end