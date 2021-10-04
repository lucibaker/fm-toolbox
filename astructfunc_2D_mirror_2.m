function [rvec,avec,avcorr_u,varargout] = astructfunc_2D_mirror_2(pos,rvec,sectors,xlim,ylim,mirror,pvar,quarters,q_idx)
% Description: ADF function of particle velocity correlation
% Equation: 
%               mean of particle-pair velocity correlation in bin (r,theta) 
% g(r,theta) = ----------------------------------------------------------------------------
%               ensemble average of particle velocity variance across all snapshots 
%                   *if one snapshot only, then particle velocity variance in the field 
% 
% Method:
% 0) Take particle positions, (if mirror - mirrors around all boundaries (up, down,
% left, right, up-left, up-right, etc.) 
% 1) Take each particle from original locations, calculate the distance of this particle to all other
% particles in the snapshot (including mirrors left and right)
% 2) Take velocity information of both particles and multiples them to find
% correlation
% 3) Bin distance information according to specified bins
% 4) Take the mean of all velocity correlations found in the bin
% 5) Normalize using ensemble average particle velcoity variance 

% Parameters:
%   pos             N x 4 matrix of particle (X,Y,U,V)
%   rvec            array of radii segments for particles
%   sectors         axial bins (24 = bins of 15deg wide)
%   xlim            vector of min,max x location
%   ylim            vector of min,max y location
%   mirror          (1) mirror or (0) no mirror
%   quarters        (1) bin data into quarters or (0) retain full 360 deg
%   q_idx           (3) u velocity (4) v velocity

%%
if ~isa(pos,'cell')
    pos2{1} = pos;
else
    pos2 = pos;
end
clear pos

%% preparing bins and bincenters
% which avcorr_2D_mirror
if quarters
    angles = linspace(0,90,sectors/4+1)';
    sector_mod = 4;
else
    angles = linspace(-180,180,sectors+1)';
    sector_mod = 1;
end

if ~mirror
    rvec = rvec(1:floor(numel(rvec)/2));
end
m = size(rvec,2);
rbin = zeros(m-1,2);
rbin(:,1) = rvec(2:end) - diff(rvec)./2;
nbin = zeros(m-1,1);

n_obs_u = zeros(m-1,sectors);
avcorr_u = zeros(m-1,sectors);
error_avcorr_u = zeros(1,2);

for t = 1:length(pos2)
    pos = pos2{t};
    
    if isempty(pos)
        continue
    end
    
    %% mirroring
    if mirror
        pos(:,1) = pos(:,1) - xlim(1);
        pos(:,2) = pos(:,2) - ylim(1);
        posmir = vertcat(pos,[pos(:,1),ylim(2)-pos(:,2)+ylim(2),pos(:,3:4)],[pos(:,1),ylim(2)-pos(:,2)-ylim(2),pos(:,3:4)]); %mirror up and down
        posmir = vertcat(posmir,[xlim(2)-pos(:,1)+xlim(2),pos(:,2),pos(:,3:4)],[xlim(2)-pos(:,1)-xlim(2),pos(:,2),pos(:,3:4)]); %mirror left and right
        posmir = vertcat(posmir,[xlim(2)-pos(:,1)+xlim(2),ylim(2)-pos(:,2)+ylim(2),pos(:,3:4)]); %mirror northeast
        posmir = vertcat(posmir,[xlim(2)-pos(:,1)-xlim(2),ylim(2)-pos(:,2)+ylim(2),pos(:,3:4)]); %mirror northwest
        posmir = vertcat(posmir,[xlim(2)-pos(:,1)+xlim(2),ylim(2)-pos(:,2)-ylim(2),pos(:,3:4)]); %mirror southeast
        posmir = vertcat(posmir,[xlim(2)-pos(:,1)-xlim(2),ylim(2)-pos(:,2)-ylim(2),pos(:,3:4)]); %mirror southwest
    else
        posmir = pos;
    end
    % figure(1);scatter(posmir(:,1),posmir(:,2),'.');daspect([1 1 1]); % check mirror (disable in parallel loop)

    %% calculate distances between particles, angles between particles, make r12 array
    np = size(posmir,1);
    r12_array = zeros(np^2,1); 
    a12_array = zeros(np^2,1); 
    u12_array = zeros(np^2,1); 
    u11_array = zeros(np^2,1);
    e12_array = zeros(np^2,1);
    lastpos = 0;
    u_error = 0.1089; %m/s estimated from a 0.5px position error
    for p=1:size(pos,1)
       pos1 = repmat(pos(p,:),size(posmir,1),1);   
       r12_array(lastpos+1:lastpos+np) = sqrt(sum((pos1 - posmir).^2,2)); 
       a12_array(lastpos+1:lastpos+np) = -atan2d(pos1(:,2)-posmir(:,2),pos1(:,1)-posmir(:,1));
       u12_array(lastpos+1:lastpos+np) = (posmir(:,q_idx)-pos1(:,q_idx)).^2; 
       e12_array(lastpos+1:lastpos+np) = u_error*(posmir(:,q_idx).*pos1(:,q_idx)); 
       u11_array(lastpos+1:lastpos+np) = pos1(:,q_idx).^2;
       lastpos = lastpos+np;   
    end

    if quarters
        a12_array(a12_array>90) = 180 - a12_array(a12_array>90);
        a12_array(a12_array<-90) = 180 + a12_array(a12_array<-90);
        a12_array(a12_array<0) = -a12_array(a12_array<0);
    end
    %% bin distances to annulus and sectors
    abin = cell(1,sectors/sector_mod);
    anbin = cell(1,sectors/sector_mod);
    for i=1:sectors/sector_mod
        r12_sector = r12_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
        u12_sector = u12_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
    %     u11_sector = u11_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
        for j=1:(m-1)
            rbin(j,2) = sum(u12_sector(bitand(r12_sector(:,1) > rvec(j),r12_sector(:,1)<rvec(j+1)),1));
            nbin(j) = sum(bitand(r12_sector(:,1) > rvec(j),r12_sector(:,1)<rvec(j+1)),1); %keyboard
    %         rbin(j,2) = rbin(j,2)/mean(u11_sector(bitand(r12_sector(:,1) > rvec(j),r12_sector(:,1)<rvec(j+1)),1)); % normalize using variance of first particles in partition
%             rbin(j,2) = rbin(j,2)/pvar; % normalize using field variance
        end
        abin{i} = rbin;
        anbin{i} = nbin;
    end

    avcorr = vertcat(abin{:});
    n_obs = vertcat(anbin{:});
    avcorr_u = avcorr_u + reshape(avcorr(:,2),[m-1,sectors]);
    n_obs_u = n_obs_u + reshape(n_obs,[m-1,sectors]);
    error_avcorr_u = error_avcorr_u + [mean(e12_array),max(e12_array)];
    
end

avec = angles(1:end-1)+0.5*mean(diff(angles));
rvec = rbin(:,1);

avcorr_u = avcorr_u./n_obs_u/pvar;
error_avcorr_u = error_avcorr_u./length(pos2)/pvar;
varargout{1} = error_avcorr_u;

end