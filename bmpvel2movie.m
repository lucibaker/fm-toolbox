% bmp series and velocity fields to movie file
% Need x, y, U1, V1, img, [centers], [radii], [pxtom], [h] already in workspace
% [...] if plotting spheres

mname = 'hgel_movie_VF10_MS8_paired_short2'; % movie file name (avi or mp4)
bname = '9991\run_f_'; % base .bmp file name

% v = VideoWriter(mname);%,'Grayscale AVI');
% v.FrameRate = 15;
% open(v)

figure(1);

N = 250; % number of frames
for i = 104%1:N
    %%
    i=i+10;
    A = imread([bname sprintf('%04d',2*i-1) '_Cam_9991_Cine1.bmp']);
    A = A(:,1.5*h.dx:h.ix-1.5*h.dx);
%     A = A(3*h.dy:(length(y)+6)*h.dy,3*h.dx:h.ix-2*h.dx);
    ax1 = subplot(121); imshow(A); colormap(ax1,'gray');
%     axis([0.004562 0.05626 0 0.149])
    hold on;
	rectangle('position',[170 320 10/1000/pxtom 5],'facecolor','white','edgecolor','none');
    annotation('textbox','position',[.261 .85 .2 .05],'string','10 mm', ...
        'fontsize',11,'fontweight','bold','edgecolor','none','color','white','margin',0);
    % ^ disable goodplot annotations
    
    ax2 = subplot(122);
    pcolor(x,y,sqrt(U1(:,:,i).^2 + V1(:,:,i).^2)); colormap(ax2,'jet'); shading flat; 
    axis equal tight; axis([0.005537 0.0443 0 0.1357])
    ca=gca; ca.CLim = [0.1, .23]; hold on; 
    c=colorbar; c.Label.String = '$u$ (m/s)';
    quiver(x,y,U1(:,:,i),V1(:,:,i),'color','k');
    if ~isempty(radii{2*i-1})
        viscircles(pxtom*centers{2*i},pxtom*radii{2*i}*1.02,'Edgecolor','b');
%         viscircles(pxtom*centers{2*i-1}+repmat([-3*pxtom*h.dx,pxtom*h.dx],length(radii{2*i-1}),1),pxtom*radii{2*i-1},'Edgecolor','b');
    end
    hold off
    
    %remove ticks, ticklabels, make axes fullscreen
    ax2.XTick = []; ax2.YTick = [];
    set(ax2,'Position',[0.4 0.05 .4 0.9],'YLim',[.002 .0858])
    set(ax1,'Position',[0 0.05 .4 0.9],'YLim',[300 800])    
%     goodplot([7 6])
    goodplot([6 5.143])
    %%

%     cdata = print('-RGBImage');
%     writeVideo(v,cdata);
end

% close(v)