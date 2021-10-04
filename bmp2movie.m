% bmp series to movie file (gif or avi)

fmt = 'gif'; % format option: 'gif' or 'avi'

mname = 'rods_hiRe_short'; % movie file name 
bname = 'images/run_a_'; % base .bmp file name
i0 = 100580; %403480; %206470; % starting frame 
N = 300;  % number of frames

if strcmp(fmt,'avi')
    v = VideoWriter(mname,'Grayscale AVI');
    v.FrameRate = 5;
    open(v)
elseif strcmp(fmt,'gif')
else
    error('format must be ''gif'' or ''avi''')
end

for i = i0:(i0 + N) 
    A = imread([bname num2str(i) '.bmp']); % '_Cam_9992_Cine1.bmp']);
    if strcmp(fmt,'avi')
        % write to video
        writeVideo(v,A);
    else
        % write to gif
        [imind,cm] = gray2ind(A,256);
        if i == i0
            imwrite(imind,cm,[mname '.gif'],'gif', 'Loopcount',inf,'DelayTime',1/10);
        else
            imwrite(imind,cm,[mname '.gif'],'gif','WriteMode','append','DelayTime',1/10);
        end
    end
    
    if ~mod(i,100)
        fprintf([num2str(i-i0) '/' num2str(N) '\n'])
    end
    
end

if strcmp(fmt,'avi')
    close(v)
end