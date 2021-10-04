function [q_lmean_bin, N_pts, t_lags] = lagrang_mean2(smtracks, q1, y_bins, y0m, Tmax, fs, pn)
% Lagrangian mean of one variable specified by q1 along particle tracks
% contained in smtracks. Options for q1 are:
%   3: u_p, 4: v_p, 8: a_xp, 9: a_yp, 10: u_f|p, 11: v_f|p, 12: theta, 13:
%   vtp, 14: atp, 15: anp
%
% Tracks are binned by y-coordinate. N_bins specifies the number of bins.
% y_bins specifies the edges of the bins on y.
%
% Tmax: maximum track length. 
%
% if pn = 1: u_lmean_bin, N_pts, t_lags contains autocorrelation and number of 
% observations for positive lags only. Size of each is Tmax x N_bins.
% if pn = -1: u_lmean_bin, N_pts, t_lags contains autocorrelation and number of 
% observations for negative lags only. Size of each is Tmax x N_bins.
% if pn = 0: u_lmean_bin, N_pts, t_lags contains autocorrelation and number of 
% observations for both positive and negative lags. Size of each is 2*Tmax-1 x N_bins.

N_bins = length(y_bins) - 1;

% lagrangian mean for each bin
%% positive lags
if pn >= 0
    
    N_pts_p = zeros(N_bins,Tmax);
    q_lmean_bin_p = zeros(N_bins,Tmax);

    for j = 1:N_bins
        for i = 1:max(smtracks(:,5))
            idx = find(smtracks(:,5)==i);
            if idx
                % times when track starts in y bin
                inbin = (smtracks(idx(1),2)-y0m >= y_bins(j) & smtracks(idx(1),2)-y0m < y_bins(j+1)); 
                if inbin
                    % accumulate track u velocity and number of samples at each time lag (aligned at t0 = 0 lag)
                    q_padded = [smtracks(idx,q1)', zeros(1, Tmax-length(idx))];
%                     q_padded2 = q_padded;
                    q_padded(isnan(q_padded)) = 0;
                    q_lmean_bin_p(j,:) = q_lmean_bin_p(j,:) + q_padded;
                    N_padded = [~isnan(q_padded(1:length(idx))), zeros(1, Tmax-length(idx))];
                    N_pts_p(j,:) = N_pts_p(j,:) + N_padded; 
                end
            end
        end   
    end
    q_lmean_bin_p = q_lmean_bin_p./N_pts_p;
    q_lmean_bin_p(isnan(q_lmean_bin_p)) = 0;
end

%% negative lags
if pn <= 0    
    
    N_pts_n = zeros(N_bins,Tmax);
    q_lmean_bin_n = zeros(N_bins,Tmax);

    for j = 1:N_bins
        for i = 1:max(smtracks(:,5))
            idx = find(smtracks(:,5)==i);
            if idx
                % times when track ends in y bin
                inbin = (smtracks(idx(end),2)-y0m >= y_bins(j) & smtracks(idx(end),2)-y0m < y_bins(j+1)); 
                if inbin
                    % accumulate track u velocity and number of samples at each time lag (aligned at t0 = 0 lag)
                    q_padded = [zeros(1, Tmax-length(idx)), smtracks(idx,q1)'];
%                     q_padded2 = q_padded;
                    q_padded(isnan(q_padded)) = 0;
                    q_lmean_bin_n(j,:) = q_lmean_bin_n(j,:) + q_padded;
                    N_padded = [zeros(1, Tmax-length(idx)), ~isnan(q_padded(1:length(idx)))];
                    N_pts_n(j,:) = N_pts_n(j,:) + N_padded; 
                end
            end
        end   
    end
    q_lmean_bin_n = q_lmean_bin_n./N_pts_n;
    q_lmean_bin_n(isnan(q_lmean_bin_n)) = 0;
end


%% outputs
if pn > 0  % positive only
    q_lmean_bin = q_lmean_bin_p;
    N_pts = N_pts_p;
    t_lags = (0:Tmax-1)/fs;

elseif pn < 0   % negative only
    q_lmean_bin = q_lmean_bin_n;
    N_pts = N_pts_n;
    t_lags = (0:Tmax-1)/fs;
    
else    % pos and neg
    q_lmean_bin = [q_lmean_bin_n, q_lmean_bin_p(:,2:end)];
    N_pts = [N_pts_n, N_pts_p(:,2:end)];
    t_lags = (-(Tmax-1):(Tmax-1))/fs;
end
    