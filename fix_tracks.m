function [tracks, tracklength] = fix_tracks(tracks, tracklength, searchrad, dt, skipframes)
% input tracks, tracklengths arrays; search radius; and maximum number of
% skipped frames to search

fprintf('\nrepairing broken tracks...')

plot_on = 0;
if plot_on
    figure;
    xlabel('x [m]'); ylabel('y [m]'); 
    c = jet(max(tracks(:,5)));
end

% fix broken tracks
for uid = 1:max(tracks(:,5))
    uid_idx = tracks(:,5) == uid;
    track1 = tracks(uid_idx,:);
    
    % endpoint of this track
%     [l_end,idx_end] = max(track1(:,6));
    l_end = track1(end,6);
    f_end = track1(end,7);
    x_end = track1(end,1);
    y_end = track1(end,2);
    u_end = track1(end,3);
    v_end = track1(end,4);
    ax_end = track1(end,8);
    ay_end = track1(end,9);

    for df = 2:skipframes
        % predicted startpoint of next track
        f_start = f_end + df;
        x_start = x_end + u_end*dt + 0.5*ax_end*dt^2;
        y_start = y_end + v_end*dt + 0.5*ay_end*dt^2;
    
        % find matching startpoints of other tracks
        match_idx = find( ~uid_idx & ... % not original track
            tracks(:,6) == 0 & ... % check for start of track
            tracks(:,7) == f_start & ... % check for frame number
            abs(tracks(:,1) - x_start) < searchrad & ... % check for x position
            abs(tracks(:,2) - y_start) < searchrad ); % check for y position
        
        if numel(match_idx) > 0

            % if more than one match, choose the one closest to the predicted startpoint
            if numel(match_idx) > 1
                [~,I] = min((tracks(match_idx,1) - x_start).^2 + (tracks(match_idx,2) - y_start).^2);
                match_idx = match_idx(I);
            end
            match_uid = tracks(match_idx,5);
            match_idx = tracks(:,5) == match_uid;
            
            % reassign uid of this track to the uid of the found match
            tracks(uid_idx,5) = match_uid;  
            
            % update matched track's lifetime counter
            tracks(match_idx,6) = ((l_end + df + 1):(l_end + df + sum(match_idx)))';

            % update tracklengths
            tracklength(match_uid) = tracklength(match_uid) + tracklength(uid);
            tracklength(uid) = nan;

            % interpolate missing values
            fixed_uid = tracks(:,5) == match_uid; % fixed track indices
            l_interp = ((l_end + 1):(l_end + df - 1))'; % frames to interpolate
            
            track_interp = nan(length(l_interp),size(tracks,2));
            track_interp(:,1:2) = interp1(tracks(fixed_uid,6),tracks(fixed_uid,1:2),l_interp); % interpolate positions
            track_interp(:,5) = match_uid*ones(size(l_interp)); % uid
            track_interp(:,6:7) = [l_interp, ((f_end+1):(f_start-1))']; % lifetime and frame counters
            track_interp(:,10:end) = interp1(tracks(fixed_uid,6),tracks(fixed_uid,10:end),l_interp); % errchk and orientations
                
            tracks = vertcat(tracks,track_interp);

            if plot_on && sum(fixed_uid) > 5
                fixed_uid = tracks(:,5) == match_uid;
                plot(tracks(fixed_uid,1),tracks(fixed_uid,2),'.','color',c(randi(length(c)),:));
                hold on; plot(tracks(find(match_idx,1),1),tracks(find(match_idx,1),2),'ro','markersize',8); hold off
                axis equal; axis([-.5 .5 -.45 .05]);
                pause; 
            end
            
            break  % when match is found, move to the next track
        end    
    end
end

tracklength(isnan(tracklength)) = [];

fprintf('found %i matches\n',uid-length(tracklength))

end