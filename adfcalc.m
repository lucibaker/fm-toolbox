function [rbin,rvec,thvec]=adfcalc(C,r_edges,theta_edges,xlims,ylims,ixlims,iylims)
% Parameters:
%   C       matrix (n,2) of n (x,y) particle locations in a single realization
%           OR cell array of (n,2) matrices containing particle locs in multiple realizations
%   rvec    row vector of radius segments for particles
%   theta   row vector of angle segments for particles (in radians) (theta =
%           angle counter-clockwise from positive x-axis)
%   xlims   vector of min,max x location of particles to consider
%   ylims   vector of min,max y location of particles to consider
%   ixlims  vector of min,max image (or domain) x boundary
%   iylims  vector of min,max image (or domain) y boundary
% 
% This function takes as input a matrix containing the co-ordinates of the
% particles, a vector of shell radii, and a vector of angles, and finds the
% ADF. 
%
% Note: theta is defined for [-pi,pi] radians.
%
% rbin(:,i) is the ADF vector, normalized by the ADF of a uniform
% distribution, for the i-th theta bin.

if ~isa(C,'cell')
    C2{1} = C;
else
    C2 = C;
end
clear C

k = size(theta_edges,2);
m = size(r_edges,2);

rbin = zeros(m-1,k-1);
n1bin = zeros(m-1,k-1);
n1 = 0;
n2 = 0;
total_np = 0;

% precompute bin edgepoints
bin_edgept = zeros(2,m-1,k-1);
for j = 1:m-1   % for each annulus
    for h = 1:k-1   % for each angle bin
        bin_edgept(:,j,h) = [r_edges(j), mean([theta_edges(h),theta_edges(h+1)])]; % inner bin edgepoint [r, theta] relative to P1
        bin_edgept(:,j,h) = [bin_edgept(1,j,h).*cos(bin_edgept(2,j,h)), bin_edgept(1,j,h).*sin(bin_edgept(2,j,h))];  % inner bin edgepoint [x, y] relative to P1
    end
end

for t = 1:length(C2)
    C = C2{t};
    total_np = total_np + size(C,1);
    
    if isempty(C)
        continue
    end
    
    % subselect set of first particles (P1)
    subselect1 = (C(:,1)>(xlims(1)) & C(:,1)<(xlims(2)) & C(:,2)<(ylims(2)) & C(:,2)>(ylims(1)));
    % subselect set of second particles (P2)
    subselect2 = (C(:,1)>(xlims(1)-r_edges(end)) & C(:,1)<(xlims(2)+r_edges(end)) & C(:,2)<(ylims(2)+r_edges(end)) & C(:,2)>(ylims(1)-r_edges(end)));
    
    x1 = C(subselect1,1);
    x2 = C(subselect2,1);
    y1 = C(subselect1,2);
    y2 = C(subselect2,2);

    clear C
    clear subselect1
    clear subselect2
    
    n1 = n1 + length(x1);
    n2 = n2 + length(x2);
    
    % precompute distances and angles
    for i = 1:length(x1)   % for all P1
        % find distance to P2
        r12 = sqrt((x1(i)-x2(:)).^2 + (y1(i)-y2(:)).^2);
        % find angle theta to P2
        th12 = atan2((y2(:)-y1(i)), (x2(:)-x1(i)));
        th12(th12 < theta_edges(1)) = th12(th12 < theta_edges(1)) + 2*pi;
        for j = 1:m-1   % for each annulus
            for h = 1:k-1   % for each angle bin 
                bin_edgept_ijh = bin_edgept(:,j,h) + [x1(i);y1(i)];  % inner bin edgepoint [x, y] in image coords
                
                % check if the edgepoint of bin [j,h] is within the image frame for this P1
                if all(bin_edgept_ijh < [ixlims(2);iylims(2)] & bin_edgept_ijh > [ixlims(1);iylims(1)])
                    % update P2 count for this bin 
                    rbin(j,h) = rbin(j,h) + sum(r12 >= r_edges(j) & r12 <= r_edges(j+1) & th12 >= theta_edges(h) & th12 <= theta_edges(h+1));
                    % increment valid P1 count for this bin
                    n1bin(j,h) = n1bin(j,h) + 1;
                end
            end
        end
    end  
end

% normalization values
total_pairs = n1bin.*(total_np-length(C2));  

for j = 1:m-1     % for each annulus
    % area or volume of bin
    area_jh = pi*((r_edges(j+1))^2-r_edges(j)^2)/(k-1);
    for h = 1:k-1   % for each angle 
        % normalize by number expected if pair distribution were random
        rbin(j,h) = rbin(j,h)/(area_jh/(diff(ixlims)*diff(iylims))*total_pairs(j,h)/length(C2)); 
    end
end

rvec = r_edges(2:end)-0.5*diff(r_edges); % radius bin centers
thvec = theta_edges(2:end)-0.5*diff(theta_edges); % angle bin centers

end