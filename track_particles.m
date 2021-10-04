%% Particle Lagrangian Tracking
% Creates tracks, smoothed tracks, and snapshots from particle position and
% velocity data.
% 
% Uses: gauss_velocity.m - Kasey Laurent
% Author : Kee Onn Fong(fongx065@umn.edu)

function [tracks,smtracks,snapshots] = track_particles(snapshots,searchrad,kernel,fs,varargin)
if ~isempty(varargin)
    u_mean = varargin{1}; % PIV predictor field
    v_mean = varargin{2};
else
    u_mean = 0;
    v_mean = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%
fprintf('finding particle tracks...')
tstart = tic;

nframes = numel(snapshots);
T = 1/fs; % T is in s
savefile=0;
%%%%%%%%%%%%%%%%%%%%%%%
%% Particle UID Assigning using knnsearch
% [X_m Y_m U_ms V_ms UID lifetime];
% tic
fA = snapshots{1};
fA(:,5) = (1:numel(fA(:,1)))'; %Unique ID assigned for first frame
fA(:,6) = 0; % lifetime of particle
snapshots{1} = fA;
maxmaxUID = 0; % overall max UID

for i = 1:nframes-1
    if isempty(snapshots{i}) || isempty(snapshots{i+1})
        continue;
    end
    fA = snapshots{i};
    fB = snapshots{i+1};
    flag = 0;
    if isempty(fB)
        disp(['frame #',num2str(i),' is empty'])
        continue;
    end
    if isempty(fA)
        count = i;
        flag = 1;
        while isempty(fA)
            fA = snapshots{count};
            count = count-1;
        end
    end
    fA = sortrows(fA,1);%fA = fA(1:100,:);
    fB = sortrows(fB,1);%fB = fB(1:100,:);
    [np,~]=size(fB);
    maxUIDfA = max([max(fA(:,5)),maxmaxUID]);
    maxmaxUID = maxUIDfA;
    minUIDfA = min(fA(:,5));
    fB(:,5) = ones(numel(fB(:,1)),1)*-1; %-1 is no ID assigned yet
    fB(:,6) = -1; % lifetime of particle 

    % Apply PIV predictor field
    fBshifted = fB;
    if numel(u_mean) > 1   
        fBshifted(:,1) = fBshifted(:,1) - T*interp2(u_mean(:,:,2),u_mean(:,:,3),u_mean(:,:,1),fBshifted(:,1),fBshifted(:,2));
        fBshifted(:,2) = fBshifted(:,2) - T*interp2(u_mean(:,:,2),u_mean(:,:,3),v_mean(:,:,1),fBshifted(:,1),fBshifted(:,2));
    end
    
    %use k-nearest neighbor search to associate particles in fA and fB
    [idx,d] = knnsearch(fA(:,1:2),fBshifted(:,1:2));
    
    idx = [idx,(1:np)',d];
    idx = sortrows(idx,3); idxuncut = numel(idx(:,1));
    idx = idx((idx(:,3)<searchrad),:);
    idx = sortrows(idx,1);

    % if a particle in frame A matches with more than 1 particle in frame B
    rowdel = [];
    for k = 1:numel(idx(:,1))
        rowno = idx(k,1);
        if sum(idx(:,1)==rowno)>1
%             disp(['duplicate found at ',num2str(rowno),' in frame ',num2str(i)]);
            rownos = find(idx(:,1)==rowno);
            rowdel = vertcat(rowdel,rownos(2:end)); % keep only the closest match
        end
    end

    rowdel = unique(rowdel);cutted = 0;
    for k=1:numel(rowdel)
        idx(rowdel(k)-cutted,:)=[];
        cutted = cutted+1;
    end
    
    %check if particle fields are time-resolved
    if flag == 1%numel(idx(:,1))<0
%         disp(['particle fields #',num2str(i),' and #',num2str(i+1),...
%            ' are NOT sequential, num(idx) = ',num2str(numel(idx(:,1)))]);
%         assign new UID for particles in fB and new lifetime count
        maxUIDfB = maxUIDfA;
        for j=1:numel(fB(:,1))
            if fB(j,5)==-1
                maxUIDfB = maxUIDfB+1;
                fB(j,5) = maxUIDfB;%+1;
                fB(j,6) = 0;
            end
        end
    else
        %disp(['particle fields of #',num2str(i),' and #',num2str(i+1),...
        %    ' are sequential']);
        %associate particles in fA with itself in fB and add lifetime count
        for j=1:numel(idx(:,1))
%             if fB(idx(j,2),1)<fA(idx(j,1),1) %% disable backwards motion
%                 continue;
%             end            
            UID = fA(idx(j,1),5);
            fB(idx(j,2),5) = UID;
            fB(idx(j,2),6) = fA(idx(j,1),6)+1;
        end
        %assign new UID for particles in fB and new lifetime count
        maxUIDfB = maxUIDfA;
        for j=1:numel(fB(:,1))
            if fB(j,5)==-1
                maxUIDfB = maxUIDfB +1;
                fB(j,5) = maxUIDfB;%+1;
                fB(j,6) = 0;
            end
        end
    end        
    
    snapshots{i+1} = fB;
end

% disp('Particle UID complete')
% toc

% %% REPAIR BROKEN TRACKS
% df = 1; % max # of frames away to look for matches
% rep_rad_x = searchrad;
% rep_rad_y = searchrad;
% 
% for i = 1:nframes-1
%     uid_f1 = snapshots{i}(:,5); % UID's in frame i
%     uid_f2 = snapshots{i+1}(:,5); % UID's in frame i+1
%     endpts = logical(~ismember(uid_f1,uid_f2)); % tracks that end in frame i
%     if any(endpts)
%         x_end = snapshots{i}(endpts,1); % x,y-coords of track endpoints
%         y_end = snapshots{i}(endpts,2); 
%         uid_end = snapshots{i}(endpts,5); 
%         lt_end = snapshots{i}(endpts,6); 
%         
%         for j = 2:(df+1)
%             startpts = logical(snapshots{i+j}(:,6) == 0);
%             if any(startpts)
%                 x_start = snapshots{i+j}(startpts,1); % x,y-coords of nearby track startpoints
%                 y_start = snapshots{i+j}(startpts,2); 
%                 uid_start = snapshots{i+j}(startpts,5); 
%                 
%                 matches = zeros(length(x_end),length(x_start));
%                 for k = 1:length(x_end)
%                     dist_x = x_start - x_end(k);
%                     dist_y = y_start - y_end(k);
%                     shift_x = j*T*interp2(u_mean(:,:,2),u_mean(:,:,3),u_mean(:,:,1),x_end(k),y_end(k));
%                     shift_y = j*T*interp2(u_mean(:,:,2),u_mean(:,:,3),v_mean(:,:,1),x_end(k),y_end(k));
%                     matches(k,:) = abs(dist_x - shift_x) < j*rep_rad_x & abs(dist_y - shift_y) < j*rep_rad_y;
%                 end
% 
%                 if any(matches)
%                     for k = 1:length(x_end)
%                         uid_start0 = uid_start(logical(matches(k,:))); % track(s) that are connected to track k ending in frame i
%                         % test if only 1 new track matches with track k, and
%                         % that new track only has one backward match which is track k
%                         if length(uid_start0) == 1 && sum(matches(:,logical(matches(k,:)))) == 1
%                             % attach tracks
%                             iter = 1;
%                             for m = (i+1):(i+j-1)
%                                 snapshots{m} = [snapshots{m}; zeros(1, size(snapshots{m},2))];
%                                 snapshots{m}(end,1) = nan; snapshots{m}(end,2) = nan; % filled centroids assigned NaN's for now
%                                 snapshots{m}(end,5) = uid_end(k);
%                                 snapshots{m}(end,6) = lt_end(k) + iter;
%                                 iter = iter + 1;
%                             end
%                             for m = (i+j):min([(i+300),nframes])
%                                 idx_snap = logical(snapshots{m}(:,5) == uid_start0);
%                                 snapshots{m}(idx_snap,5) = uid_end(k);
%                                 snapshots{m}(idx_snap,6) = snapshots{m}(idx_snap,6) + lt_end(k) + 1 + j;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end 
% end

        

%% New Particle Velocity
% recalculate velocity for connected particles

% tic
snapminUID = zeros(nframes,1);
snapmaxUID = zeros(nframes,1);
for i=1:nframes-1
    if isempty(snapshots{i}) || isempty(snapshots{i+1})
        snapminUID(i) = max(snapminUID);
        snapmaxUID(i) = max(snapmaxUID);
        continue;
    end
    fA = snapshots{i};
    fB = snapshots{i+1};
    for fbx=1:numel(fB(:,1))
       if fB(fbx,6)>0
          UID = fB(fbx,5);
          fax = find(fA(:,5)==UID);
          U_ms = (fB(fbx,1)-fA(fax,1))/T;
          V_ms = (fB(fbx,2)-fA(fax,2))/T;
          fA(fax,3) = U_ms;
          fA(fax,4) = V_ms;
       end
    end  
    vels = fA(:,3:4);
    vels(vels==0) = nan; % remove zero velocities on track ends
    fA(:,3:4) = vels;
    snapshots{i} = fA;
    snapminUID(i) = min(snapshots{i}(:,5));
    snapmaxUID(i) = max(snapshots{i}(:,5));
end
if isempty(snapshots{nframes})
    snapminUID(nframes) = max(snapminUID);
    snapmaxUID(nframes) = max(snapmaxUID);
else
    snapminUID(nframes) = min(snapshots{nframes}(:,5));
    snapmaxUID(nframes) = max(snapshots{nframes}(:,5));
end
% plot(snapminUID); hold on; plot(snapmaxUID);
% disp('Particle Velocity Correction Complete')
% toc
%% Forming Particle Tracks
% Previously all particle locations are sorted by frame number, 
% now they are going to be sorted by UID - particle tracks
% [X_m Y_m U_ms V_ms UID lifetime frame#];

ntracks = maxUIDfB;
tracks = cell(1,ntracks);
notempty = ones(1,ntracks);

% tic
for i=1:ntracks
    snapmin = find(snapmaxUID>=i,1,'first');
    snapmax = find(snapminUID<=i,1,'last');
    if snapmax < snapmin    % track is empty (only one particle)
        notempty(i) = 0;
    else
        for j=snapmin:snapmax % frame number
            if ~isempty(snapshots{j})
                index = find(snapshots{j}(:,5)==i); % locate UID in frame
                if ~isempty(index)  
                    if size(snapshots{j},2) > 6
                        tracks{i} = vertcat(tracks{i},[snapshots{j}(index,1:6),j,snapshots{j}(index,7:8)]);
                    else
                        tracks{i} = vertcat(tracks{i},[snapshots{j}(index,:),j]);
                    end
                end
            end
        end
        if ~isempty(tracks{i})
            Udot_ms2 = gradient(tracks{i}(:,3))/T;
            Vdot_ms2 = gradient(tracks{i}(:,4))/T;
        else
            Udot_ms2 = [];
            Vdot_ms2 = [];
        end
        if size(tracks{i},2) > 7
            tracks{i} = [tracks{i}(:,1:7) Udot_ms2 Vdot_ms2 tracks{i}(:,8:9)];
        else
            tracks{i} = [tracks{i} Udot_ms2 Vdot_ms2];
        end
        tracks{i}(:,7) = tracks{i}(:,7) - 2;
    end
end
tracks = tracks(logical(notempty));  % remove empty tracks
ntracks = length(tracks);

% disp('Forming Particle Tracks Complete');
% toc    

%% Particle Velocity Smoothing and Acceleration Calculation
% Acceleration = (U2-U1)/T in three directions 
% Vdot_ms2 is not valid due to resolution
% only calculated for particles with long lifetime to be meaningful
% gauss_velocity.m is used here
% [X_mm Y_mm U_ms V_ms UID lifetime frame# Udot_ms2 Vdot_ms2];

smtracks = cell(1,ntracks);

if kernel > 1
    mintracks = kernel+3; %minimum length of tracks for smoothing (at least three points)
    % tic
    for i=1:ntracks
        if length(tracks{i}(:,1))<mintracks
            continue;
        else
            xx = tracks{i}(:,1);
            yy = tracks{i}(:,2);
            [x2, range] = gauss_position(xx',kernel);
            [y2, ~] = gauss_position(yy',kernel);
            [u2, ~] = gauss_velocity(xx',kernel,T);  % combo smoothing + differentiating filter
            [v2, ~] = gauss_velocity(yy',kernel,T);

            UID2 = ones(1,length(x2))*i;
            if ~isempty(x2)
                [Udot_ms2,~] = gauss_accel(xx',kernel,T);
                [Vdot_ms2,~] = gauss_accel(yy',kernel,T);       
                if size(tracks{i},2) > 9
                    smtracks{i} = [x2',y2',u2',v2',UID2',range', ...
                    range'+tracks{i}(1,7),Udot_ms2',Vdot_ms2',tracks{i}(range,10:11)];
                else
                    smtracks{i} = [x2',y2',u2',v2',UID2',range', ...
                    range'+tracks{i}(1,7),Udot_ms2',Vdot_ms2'];
                end
            end
        end
    end
else
    mintracks = 1;
    for i=1:ntracks
        if length(tracks{i}(:,1))<mintracks
            continue;
        else
            smtracks{i} = tracks{i};
        end
    end
end

% disp('Particle Smoothing and Acceleration Calc Complete')
% toc

%% save UID file
if savefile
    save([load_str,'_UID',num2str(kernel),'.mat'],'snapshots','tracks','smtracks','kernel');
    disp('Track File Saved')
end

fprintf('done. ')
toc(tstart)

end