function img_to_gif(img_fn,gif_fn,delay,range)
% convert .img series to .gif
% img_fn: .img filename (including extension)
% gif_fn: .gif filename (including extension)
% delay: interval between frames [s] (1/framerate)
% range: (optional) 2-element vector [startframe, endframe] starting from 0

addpath('Z:\luci\MATLAB\pivio-matlab','Z:\luci\MATLAB\PIV-PTV-toolbox')
img = imgRead(img_fn);

if nargin < 4
    range = [0,img.it-1];
end

figure;
% set(gcf, 'position', [680   608   560   370]);
set(gcf, 'position', [683   600   560   342]);

for i = range(1):2:range(2)
    % image frame
    pcolor_img(flipud(imgGetFrame(img,i)));
    axis([0 1182 0 714]);
    set(gca,'YTick',[]); set(gca,'XTick',[]);
    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                     %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
    
    pause(1/100)

    % write to gif
    cdata = print('-RGBImage','-r300');
    frame.cdata = cdata; frame.colormap = [];
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == range(1)
        imwrite(imind,cm,gif_fn,'gif','DelayTime',delay,'Loopcount',inf);
    else
        imwrite(imind,cm,gif_fn,'gif','DelayTime',delay,'WriteMode','append');
    end
end