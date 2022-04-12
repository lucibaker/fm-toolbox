function imgset_to_gif(imgset,delay,range)
% convert image series to .gif
% imgset: image set attributes variable returned by dir(<image set directory>)
% gif_fn: .gif filename (including extension)
% delay: interval between frames [s] (1/framerate)
% range: (optional) 2-element vector [startframe, endframe] starting from 0

if nargin < 4
    range = [1,length(imgset)];
end

figure;
set(gcf, 'position', [683   600   560   342]);

for i = range(1):2:range(2)
    % image frame
    pcolor_img(flipud(imgGetFrame(img,i)));
    set(gca,'YTick',[]); set(gca,'XTick',[]);
    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                     %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
    
    pause(1/100)

    % write to gif
    fig_to_gif(gif_fn,delay)
end