% velocity fields to movie file
% Need x, y, U1, V1, img, [centers], [radii], [pxtom], [h] already in workspace
% [...] if plotting spheres

mname = 'hgel_movie_VF10_MS8_vel_short.avi'; % movie file name (avi or mp4)

v = VideoWriter(mname);%,'Grayscale AVI');
v.FrameRate = 20;
open(v)

figure(3);

N = 1000;%h.nt; % number of frames
for i = 1:N
    pcolor(x,y,sqrt(U1(:,:,i).^2 + V1(:,:,i).^2)); colormap('jet'); shading flat; 
    axis equal; axis([0 0.0319 0 0.0943])
    ca=gca; ca.CLim = [0.1, .24]; hold on; 
%     quiver(x,y,U1(:,:,i),V1(:,:,i),'color','k');
    if ~isempty(radii{2*i-1})
        viscircles(pxtom*centers{2*i-1}+repmat([-3*pxtom*h.dx,pxtom*h.dx],length(radii{2*i-1}),1),pxtom*radii{2*i-1},'Edgecolor','b');
    end
    hold off
    F = getframe;
%     cdata = print('-RGBImage');
    writeVideo(v,F.cdata);
end

close(v)