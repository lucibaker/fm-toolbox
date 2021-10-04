function [r, N_pts, t_lags] = lagrang_xcorr(smtracks, q1, q2, y_bins, y0m, Tmax, fs, pn)
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
r_p_num = zeros(N_bins,Tmax);
r_p_den1 = zeros(N_bins,1);
r_p_den2 = zeros(N_bins,Tmax);
N_pts_p = zeros(N_bins,Tmax);

for j = 1:N_bins
    for i = 1:max(smtracks(:,5))
        idx = find(smtracks(:,5)==i);
        if length(idx) > 2
            % times when track is in bin
            inbin = find( (smtracks(idx,2)-y0m >= y_bins(j) & smtracks(idx,2)-y0m < y_bins(j+1)) ); 
            if ~isempty(inbin) % if track i is in bin j
                for k = 1:length(inbin)
                    % calculate autocorrelation of [u of track i] - [lagrangian mean velocity of bin j]
                    ufluct = smtracks(idx(inbin(k):end),q1)' - u_lmean_bin_p(j,1:length(idx(inbin(k):end)));
                    vfluct = smtracks(idx(inbin(k):end),q2)' - v_lmean_bin_p(j,1:length(idx(inbin(k):end)));
                    r_num = ufluct(1)*vfluct;  % numerator of correlation coefficient
                    r_den1 = ufluct(1);     % first term of denominator
                    r_den2 = vfluct;        % second term of denominator
                    % accumulate autocorr and number of samples at each time lag
                    r_p_num(j,:) = r_p_num(j,:) + [r_num, zeros(1, Tmax-length(r_num))];
                    r_p_den1(j) = r_p_den1(j) + r_den1^2;
                    r_p_den2(j,:) = r_p_den2(j,:) + [(r_den2.^2), zeros(1, Tmax-length(r_num))];
                    N_padded = [ones(1, length(r_num)), zeros(1, Tmax-length(r_num))];
                    N_pts_p(j,:) = N_pts_p(j,:) + N_padded;
                end
            end
        end
    end
    r_p_num(j,:) = r_p_num(j,:)./N_pts_p(j,:); % take the mean of the numerator and denominators (dividing by # observations at each lag)
    r_p_den1(j) = r_p_den1(j)/N_pts_p(j,1);
    r_p_den2(j,:) = r_p_den2(j,:)./N_pts_p(j,:);
    
    r_p(j,:) = r_p_num(j,:)./(r_p_den1(j)*r_p_den2(j,:)).^(1/2);   % normalize to obtain the correlation coefficient
end

if pn
    
    %% negative lags
    % lagrangian mean velocity for each bin
    [u_lmean_bin_n, ~, ~] = lagrang_mean(smtracks, q1, y_bins, y0m, Tmax, fs, -1);
    [v_lmean_bin_n, ~, ~] = lagrang_mean(smtracks, q2, y_bins, y0m, Tmax, fs, -1);

    % lagrangian autocorrelation for each bin
    r_n = zeros(N_bins,Tmax);
    r_n_num = zeros(N_bins,Tmax);
    r_n_den1 = zeros(N_bins,1);
    r_n_den2 = zeros(N_bins,Tmax);
    N_pts_n = zeros(N_bins,Tmax);

    for j = 1:N_bins
        for i = 1:max(smtracks(:,5))
            idx = find(smtracks(:,5)==i);
            if length(idx) > 2
                % times when track ends in y bin
                inbin = find( (smtracks(idx,2)-y0m >= y_bins(j) & smtracks(idx,2)-y0m < y_bins(j+1)) ); 
                if ~isempty(inbin)
                    for k = 1:length(inbin)
                        % calculate autocorrelation of [u of track i] - [lagrangian mean velocity of bin j]
                        ufluct = smtracks(idx(inbin(k):end),q1)' - u_lmean_bin_n(j,end-length(idx(inbin(k):end))+1:end);
                        vfluct = smtracks(idx(inbin(k):end),q2)' - v_lmean_bin_n(j,end-length(idx(inbin(k):end))+1:end); 
                        r_num = ufluct(end)*vfluct(end:-1:1);  % numerator of correlation coefficient
                        r_den1 = ufluct(end);     % first term of denominator
                        r_den2 = vfluct(end:-1:1);        % second term of denominator
                        % accumulate autocorr and number of samples at each time lag
                        r_n_num(j,:) = r_n_num(j,:) + [r_num, zeros(1, Tmax-length(r_num))];
                        r_n_den1(j) = r_n_den1(j) + r_den1^2;
                        r_n_den2(j,:) = r_n_den2(j,:) + [(r_den2.^2), zeros(1, Tmax-length(r_num))];
                        N_padded = [ones(1, length(r_num)), zeros(1, Tmax-length(r_num))];
                        N_pts_n(j,:) = N_pts_n(j,:) + N_padded;
                    end
                end
            end
        end
        r_n_num(j,:) = r_n_num(j,:)./N_pts_n(j,:); % take the mean of the numerator and denominators (dividing by # observations at each lag)
        r_n_den1(j) = r_n_den1(j)/N_pts_n(j,1);
        r_n_den2(j,:) = r_n_den2(j,:)./N_pts_n(j,:);

        r_n(j,:) = r_n_num(j,:)./(r_n_den1(j)*r_n_den2(j,:)).^(1/2);   % normalize to obtain the correlation coefficient
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
