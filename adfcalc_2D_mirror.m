function [rvec,avec,adf] = adfcalc_2D_mirror(pos,rvec,sectors,xlim,ylim,mirror,quarters)
% Description: ADF function of particle distribution using pair formula
% Equation: 
%               num pair found at sector (r,theta)            total volume
% g(r,theta) = -------------------------------------- x ------------------------------------
%               num of all possible pairs                search volume at sector (r,theta)         
% Ref: Salazar JFM 2008 https://doi.org/10.1017/S0022112008000372
% Ref: Jong IJMF 2010 https://doi.org/10.1016/j.ijmultiphaseflow.2009.11.008
% Method:
% 0) Take particle positions, mirror them around all boundaries (up, down,
% left, right, up-left, up-right, etc.) 
% 1) Take each particle from original locations, calculate the distance of this particle to all other
% particles in the snapshot (including mirrors left and right)
% 2) Consolidate all particle distances calculated into a single matrix
% 3) Bin distance information according to specified bins
% 4) Normalize using equation above (depending on dimensionality) 

% Parameters:
%   pos             2 x N matrix of particle (X,Y)
%   rvec            array of radii segments for particles
%   sectors         axial bins (24 = bins of 15deg wide)
%   xlim            vector of min,max x location
%   ylim            vector of min,max y location
%   mirror          (1) mirror or (0) no mirror
%   quarters        (1) bin data into quarters or (0) retain full 360 deg

%% preparing bins and bincenters
% which adfcalc_2D_mirror
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
r12_array = []; a12_array = [];

%% mirroring particles if required
if mirror
    posmir = vertcat(pos,[pos(:,1),ylim-pos(:,2)+ylim],[pos(:,1),ylim-pos(:,2)-ylim]); %mirror up and down
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,pos(:,2)],[xlim-pos(:,1)-xlim,pos(:,2)]); %mirror left and right
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,ylim-pos(:,2)+ylim]); %mirror northeast
    posmir = vertcat(posmir,[xlim-pos(:,1)-xlim,ylim-pos(:,2)+ylim]); %mirror northwest
    posmir = vertcat(posmir,[xlim-pos(:,1)+xlim,ylim-pos(:,2)-ylim]); %mirror southeast
    posmir = vertcat(posmir,[xlim-pos(:,1)-xlim,ylim-pos(:,2)-ylim]); %mirror southwest
else
    posmir = pos;
end
% % % figure(111);scatter(posmir(:,1),posmir(:,2),'.');daspect([1 1 1]); % check mirror (disable in parallel loop)

%% calculate distances between particles, angles between particles, make r12 array
for p=1:size(pos,1)
   pos1 = repmat(pos(p,:),size(posmir,1),1);
   r12_array = vertcat(r12_array,sqrt(sum((pos1 - posmir).^2,2))); 
   a12_array = vertcat(a12_array,-atan2d(pos1(:,1)-posmir(:,1),pos1(:,2)-posmir(:,2)));   
end

if quarters
    a12_array(a12_array>90) = 180 - a12_array(a12_array>90);
    a12_array(a12_array<-90) = 180 + a12_array(a12_array<-90);
    a12_array(a12_array<0) = -a12_array(a12_array<0);
    exactly45 = sum(a12_array==45);
    if exactly45>0; disp(['warning: exactly 45 deg not counted, N = ',num2str(exactly45)]); end
end
%% bin distances to annulus and sectors
% avg_pd = size(pos,1)/(xlim*ylim); % average particle concentrations
rbin = zeros(m-1,2); 
rbin(:,1) = rvec(2:end) - diff(rvec)./2; %radial bin centers

if quarters % for quartering, 'less than or equal to' required to be equal with RDF 
    for i=1:(sectors/sector_mod)
        r12_sector = r12_array(bitand(a12_array(:,1) >= angles(i),a12_array(:,1) <= angles(i+1)),:); 
        for j=1:(m-1)
            rbin(j,2) = sum(bitand(r12_sector(:,1) >= rvec(j),r12_sector(:,1)<= rvec(j+1)));
    %         rbin(j,2) = rbin(j,2)./(pi*(rvec(j+1).^2-rvec(j).^2)*avg_pd*size(pos,1)/(sectors/sector_mod)); %old code
            search_vol = pi*(rvec(j+1).^2-rvec(j).^2)/(sectors/sector_mod);
            total_vol = xlim*ylim;
            total_pairs = (size(pos,1)*(size(pos,1)-1));
            rbin(j,2) = (rbin(j,2)/total_pairs)*(total_vol/search_vol); % edit 
        end
        abin{i} = rbin;
    end
else % no quartering, 'less than' only
    for i=1:(sectors/sector_mod)
        r12_sector = r12_array(bitand(a12_array(:,1) >= angles(i),a12_array(:,1) < angles(i+1)),:); 
        for j=1:(m-1)
            rbin(j,2) = sum(bitand(r12_sector(:,1) >= rvec(j),r12_sector(:,1)<rvec(j+1)));
            search_vol = pi*(rvec(j+1).^2-rvec(j).^2)/(sectors/sector_mod);
            total_vol = xlim*ylim;
            total_pairs = (size(pos,1)*(size(pos,1)-1));
            rbin(j,2) = (rbin(j,2)/total_pairs)*(total_vol/search_vol); 
        end
        abin{i} = rbin;
    end
end

rvec = rbin(:,1);
avec = angles(1:end-1)+0.5*mean(diff(angles));
adf = vertcat(abin{:}); adf = reshape(adf(:,2),[numel(rvec),numel(avec)]);

clear r12_array avg_pd m a12_array posmir
end