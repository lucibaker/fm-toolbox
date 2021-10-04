function [r, N_pts, t_lags] = lagrang_xcorr2(smtracks, q1, q2, y_bins, y0m, Tmax, fs, pn)
% Lagrangian autocorrelation of one or two variables along particle tracks
% contained in smtracks. Variables specified by q1, q2; can be the same or 
% different. Options for q1, q2 are:
%   3: u_p, 4: v_p, 8: a_xp, 9: a_yp, 10: u_f|p, 11: v_f|p, 12: theta, 13:
%   vtp, 14: atp, 15: anp
%
% Tracks are binned by y-coordinate. N_bins specifies the number of bins.
% y_bins specifies the edges of the bins on y.
%
% Tmax: maximum track length. fs: sampling frequency.
%
% if pn = 0: r_up, N_pts, t_lags contains autocorrelation and number of 
% observations for positive lags only. Size of each is Tmax x N_bins.
% if pn = 1: r_up, N_pts, t_lags contains autocorrelation and number of 
% observations for both positive and negative lags. Size of each is 2*Tmax-1 x N_bins.

N_bins = length(y_bins) - 1;

%% positive lags
% lagrangian mean velocity for each bin
[u_lmean_bin_p, ~, ~] = lagrang_mean(smtracks, q1, y_bins, y0m, Tmax, fs, 1);
[v_lmean_bin_p, ~, ~] = lagrang_mean(smtracks, q2, y_bins, y0m, Tmax, fs, 1);

% lagrangian autocorrelation for each bin
r_p = zeros(N_bins,Tmax);
N_pts_p = zeros(N_bins,Tmax);

for j = 1:N_bins
    for i = 1:max(smtracks(:,5))
        idx = find(smtracks(:,5)==i);
        if length(idx) > 2
            % times when track starts in y bin
            inbin = (smtracks(idx(1),2)-y0m >= y_bins(j) & smtracks(idx(1),2)-y0m < y_bins(j+1)); 
            if inbin
                % calculate autocorrelation of [u of track i] - [lagrangian mean velocity of bin j]
                ufluct = smtracks(idx,q1) - u_lmean_bin_p(j,1:length(idx))';
                vfluct = smtracks(idx,q2) - v_lmean_bin_p(j,1:length(idx))';
                r = xcorr(ufluct,vfluct,'unbiased'); 
                r = r((end+1)/2:end); 
                % accumulate autocorr and number of samples at each time lag
                r_p(j,:) = r_p(j,:) + [r', zeros(1, Tmax-length(r))];
                N_padded = [ones(1, length(r)), zeros(1, Tmax-length(r))];
                N_pts_p(j,:) = N_pts_p(j,:) + N_padded;
            end
        end
    end
    r_p(j,:) = r_p(j,:)./N_pts_p(j,:); % normalize by # observations at each lag
    r_p(j,:) = r_p(j,:)/r_p(j,1);   % normalize by zero lag value
end

if pn
    
    %% negative lags
    % lagrangian mean velocity for each bin
    [u_lmean_bin_n, ~, ~] = lagrang_mean(smtracks, q1, y_bins, y0m, Tmax, fs, -1);
    [v_lmean_bin_n, ~, ~] = lagrang_mean(smtracks, q2, y_bins, y0m, Tmax, fs, -1);

    % lagrangian autocorrelation for each bin
    r_n = zeros(N_bins,Tmax);
    N_pts_n = zeros(N_bins,Tmax);

    for j = 1:N_bins
        for i = 1:max(smtracks(:,5))
            idx = find(smtracks(:,5)==i);
            if length(idx) > 2
                % times when track ends in y bin
                inbin = (smtracks(idx(end),2)-y0m >= y_bins(j) & smtracks(idx(end),2)-y0m < y_bins(j+1)); 
                if inbin
                    % calculate autocorrelation of [u of track i] - [lagrangian mean velocity of bin j]
                    ufluct = smtracks(idx,q1) - u_lmean_bin_n(j,end-length(idx)+1:end)';
                    vfluct = smtracks(idx,q2) - v_lmean_bin_n(j,end-length(idx)+1:end)';
                    r = xcorr(ufluct(end:-1:1),vfluct(end:-1:1),'unbiased'); 
                    r = r((end+1)/2:end); 
                    % accumulate autcorr and number of samples at each time lag
                    r_n(j,:) = r_n(j,:) + [r', zeros(1, Tmax-length(r))];
                    N_padded = [ones(1, length(r)), zeros(1, Tmax-length(r))];
                    N_pts_n(j,:) = N_pts_n(j,:) + N_padded;
                end
            end
        end
        r_n(j,:) = r_n(j,:)./N_pts_n(j,:); % normalize by # observations at each lag
        r_n(j,:) = r_n(j,:)/r_n(j,1);   % normalize by zero lag value
    end

    % combine pos and neg lags
    r = [r_n(:,end:-1:1), r_p(:,2:end)];
    N_pts = [N_pts_n(:,end:-1:1), N_pts_p(:,2:end)];
    t_lags = (-(Tmax-1):(Tmax-1))/fs;

else
    r = r_p;
    N_pts = N_pts_p;
    t_lags = (0:Tmax-1)/fs;
end


end
