function fig_to_gif(gif_fn,delay)
% write current figure as a frame of a gif
% gif_fn: .gif filename (including extension)
% delay: interval between frames [s] (1/framerate)
 
cdata = print('-RGBImage','-r450');
frame.cdata = cdata; frame.colormap = [];
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if ~exist(gif_fn,'file')
    imwrite(imind,cm,fn_gif,'gif','Loopcount',inf,'DelayTime',delay);
else
    imwrite(imind,cm,fn_gif,'gif','WriteMode','append','DelayTime',delay);
end