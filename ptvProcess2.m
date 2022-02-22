function [smtracks,smtracklength] = ptvProcess2(centers,kernel,pxtom,fs,searchrad,varargin)
% Track particles using kNN search. No PIV predictor
% (using dir, not pivio)
%
% Inputs:
% h: PIV header
% centers: cell array containing the coordinates of the particle centroids,
%   one cell per frame [px]
% kernel: smoothing kernel size
% pxtom: meters/pixel
% fs: sampling frequency of particle centroids
% searchrad: search radius [m]
% (optional) angles:  cell array containing the orientation information of 
%   the particle centroids, one cell per frame
%
% Outputs:
% smtracks: smoothed particle tracks in the format [X_m Y_m U_ms V_ms UID lifetime frameno Ax_ms2 Ay_ms2 (pitch_rad) (minor_axis_m)]
% smtracklength: length of each track (# frames)

if ~isempty(varargin)
    angles = varargin{1};
else
    angles = cell(size(centers));
end

img_nt = length(centers);

snapshots = cell(img_nt,1); 
for i = 1:img_nt 
    snapshots{i} = [centers{i}*pxtom zeros(size(centers{i},1),4) angles{i}];
end

avar_k = zeros(length(kernel),1);
for n = 1:length(kernel)
    [tracks,smtracks,snapshots] = track_particles(snapshots,searchrad,kernel(n),fs);

    %% View UID Particle Track and Lagrangian Stats
    filled = smtracks(~cellfun('isempty',smtracks))';
    if isempty(filled)
        error('no tracks found')
    end
    for i=1:length(filled)
        smUID(i,1) = mode(filled{i}(:,5));
        smUID(i,2) = length(filled{i}(:,5));
    end
    smUID = sortrows(smUID,2);
    smUID = flipud(smUID);
    
    tracklength = zeros(length(tracks),1);
    fprintf(['K = ' num2str(kernel(n)) ', Ntracks = ' num2str(length(filled)) '\n'])
    for i = 1:length(smtracks)
        tracklength(i) = size(smtracks{i},1);
    end   
    filled_tracklength = tracklength(~cellfun('isempty',smtracks));
    
    filled_array = zeros(0,size(filled{1},2));
    for i = 1:length(filled)
        filled_array = [filled_array; filled{i}];
    end

    %% x-acceleration variance
    avar_k(n) = nanstd(sqrt(filled_array(:,8).^2+ filled_array(:,9).^2))^2; % 
end

%% test kernel size with x-acceleration variance
if length(kernel) > 1
    
    figure; semilogy(kernel,avar_k,'b+','linewidth',1.5); hold on
    xlabel('$t_k^+$'); ylabel('var($a_{x,p}$)'); %title('Variance of d^2X_p/dt^2');
%     n0 = find(kernel==17); % index of minimum kernel size for good exponential fit
%     P = polyfit(kernel(n0:end)',log(avar_k(n0:end)),1); 
%     semilogy(kernel, exp(P(1)*kernel + P(2)), 'k-','linewidth',1);
%     semilogy(kernel(n0),avar_k(n0),'ro','linewidth',1.5,'markersize',8);    
    keyboard
else
    smtracks = filled_array;
    smtracklength = filled_tracklength;
end