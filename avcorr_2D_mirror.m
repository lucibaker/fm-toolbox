function [rvec,avec,avcorr_u,varargout] = avcorr_2D_mirror(pos,rvec,sectors,xlim,ylim,mirror,pvar,quarters)
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
%   pos             4 x N matrix of particle (X,Y,U,V)
%   rvec            array of radii segments for particles
%   sectors         axial bins (24 = bins of 15deg wide)
%   xlim            vector of min,max x location
%   ylim            vector of min,max y location
%   mirror          (1) mirror or (0) no mirror
%   quarters        (1) bin data into quarters or (0) retain full 360 deg

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

%% mirroring
if mirror
    posmir = vertcat(pos,[pos(:,1),ylim-pos(:,2)+ylim,pos(:,3:4)],[pos(:,1),ylim-pos(:,2)-ylim,pos(:,3:4)]); %mirror up and down
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,pos(:,2),pos(:,3:4)],[xlim-pos(:,1)-xlim,pos(:,2),pos(:,3:4)]); %mirror left and right
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,ylim-pos(:,2)+ylim,pos(:,3:4)]); %mirror northeast
    posmir = vertcat(posmir,[xlim-pos(:,1)-xlim,ylim-pos(:,2)+ylim,pos(:,3:4)]); %mirror northwest
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,ylim-pos(:,2)-ylim,pos(:,3:4)]); %mirror southeast
    posmir = vertcat(posmir,[xlim-pos(:,1)-xlim,ylim-pos(:,2)-ylim,pos(:,3:4)]); %mirror southwest
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
   a12_array(lastpos+1:lastpos+np) = -atan2d(pos1(:,1)-posmir(:,1),pos1(:,2)-posmir(:,2));
   u12_array(lastpos+1:lastpos+np) = (posmir(:,4).*pos1(:,4)); 
   e12_array(lastpos+1:lastpos+np) = u_error*(posmir(:,4).*pos1(:,4)); 
   u11_array(lastpos+1:lastpos+np) = pos1(:,4).^2;
   lastpos = lastpos+np;   
end

if quarters
    a12_array(a12_array>90) = 180 - a12_array(a12_array>90);
    a12_array(a12_array<-90) = 180 + a12_array(a12_array<-90);
    a12_array(a12_array<0) = -a12_array(a12_array<0);
end
%% bin distances to annulus and sectors
abin = cell(1,sectors/sector_mod);
for i=1:sectors/sector_mod
    r12_sector = r12_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
    u12_sector = u12_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
%     u11_sector = u11_array(bitand(a12_array(:,1) > angles(i),a12_array(:,1) < angles(i+1)),:); 
    for j=1:(m-1)
        rbin(j,2) = mean(u12_sector(bitand(r12_sector(:,1) > rvec(j),r12_sector(:,1)<rvec(j+1)),1));
%         rbin(j,2) = rbin(j,2)/mean(u11_sector(bitand(r12_sector(:,1) > rvec(j),r12_sector(:,1)<rvec(j+1)),1)); % normalize using variance of first particles in partition
        rbin(j,2) = rbin(j,2)/pvar; % normalize using field variance
    end
    abin{i} = rbin;
end

rvec = rbin(:,1);
avec = angles(1:end-1)+0.5*mean(diff(angles));
avcorr = vertcat(abin{:}); 
avcorr_u = reshape(avcorr(:,2),[numel(rvec),numel(avec)]);
varargout{1} = [mean(e12_array)/pvar,max(e12_array)/pvar];

end