% bmp series to movie file (gif or avi)

fmt = 'avi'; % format option: 'gif' or 'avi'

mname = 'test_nurdles_raw'; % movie file name 
bname = 'Basler_acA2040-90um__23703425__20220113_100116310_'; % base .bmp file name
i0 = 1200;  % starting frame 
N = 100;  % number of frames

if strcmp(fmt,'avi')
    v = VideoWriter(mname,'Grayscale AVI');
    v.FrameRate = 10;
    open(v)
elseif strcmp(fmt,'gif')
else
    error('format must be ''gif'' or ''avi''')
end

for i = i0:(i0 + N) 
    A = imread(sprintf('%s%04i.tiff',bname, i)); 
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