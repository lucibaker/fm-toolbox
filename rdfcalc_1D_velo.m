function [rbin]=rdfcalc_1D_velo(C,V,rvecarray,xlims,pvar,varargin)
% Description: 1-D RDF function and velocity correlation using pair formula
% Equation: 
%         num pair found at dist r            total volume
% g(r) = -------------------------- x ----------------------------
%               search volume           num of all possible pairs
% Ref: Salazar JFM 2008 https://doi.org/10.1017/S0022112008000372
% Ref: Jong IJMF 2010 https://doi.org/10.1016/j.ijmultiphaseflow.2009.11.008
% Method:
% 0) Take particle positions, mirror them left and right
% 1) Take each particle from original locations, calculate the distance of this particle to all other
% particles in the snapshot (including mirrors left and right)
% 2) Consolidate all particle distances calculated into a single matrix
% 3) Bin distance information according to specified bins
% 4) Normalize using equation above (depending on dimensionality) - need to
% divide by 2 for mirroring? (tested with random generated particles)

% Input:
%   C           vector of particle x locations
%   V           vector of particle velocities associated with C
%   rvecarray   array of radii bin edges for particles
%   xlim        vector of min,max x location
%   pvar        value to normalize the velocity correlations (I use the variance of V(:))
%   varargin(1) mirroring switch (1 for on, 0 for off)

% Output:
%   rbin(:,1)     vector r to plot data (bincenters of rvecarray)
%   rbin(:,2)     RDF result 
%   rbin(:,3)     R_VV result
 
if ~isempty(varargin)
    mirror = varargin{1}; m = size(rvecarray,2);
    rbin = zeros(m-1,2);
    rbin(:,1)=rvecarray(2:end)-diff(rvecarray)/2;
else 
    mirror = 0;
    m = ceil(size(rvecarray,2)/2);
    rbin = zeros(m-1,2);
    rbin(:,1)=rvecarray(2:m)-diff(rvecarray(1:m))/2;
end

r12_array = []; u12_array = [];

if mirror
    disp('mirroring')
    for i=1:numel(C) % for each particle in frame
        C_mirror1 = (xlims(2)-C(:))-xlims(2); %mirror of particles on left
        C_mirror2 = (xlims(2)-C(:))+xlims(2); %mirror of particles on right
        r12 = abs(C(i)-C(:)); % calc distance to other particles including itself
        r12m1 = abs(C(i)-C_mirror1(:)); % calc distance to other particles in left mirror
        r12m2 = abs(C(i)-C_mirror2(:)); % calc distance to other particles in right mirror
        r12_array = vertcat(r12_array,r12,r12m1,r12m2); % store pair distance information
        u12 = V(i)*V(:); % calc distance to other particles including itself
        u12m1 = abs(C(i)-C_mirror1(:)); % calc distance to other particles in left mirror
        u12m2 = abs(C(i)-C_mirror2(:)); % calc distance to other particles in right mirror
        u12_array = vertcat(u12_array,u12,u12m1,u12m2); % store pair distance information
    end
else
    disp('no mirroring')
    for i=1:numel(C) % for each particle in frame
        r12 = abs(C(i)-C(:)); % calc distance to other particles including itself
        u12 = V(i)*V(:); %calc velocity correlation of particle pairs including itself
        r12_array = vertcat(r12_array,r12); % store pair distance information
        u12_array = vertcat(u12_array,u12); % store pair velocorr information
    end
end

for j = 1:m-1   % count for each annulus how many particle pairs are within the limits
    rbinall(j) = sum(bitand(r12_array>rvecarray(j),r12_array<rvecarray(j+1))); 
    ubinall(j) = nanmean(u12_array((bitand(r12_array>rvecarray(j),r12_array<rvecarray(j+1))))); 
end

search_volume = 2*mean(diff(rvecarray));
total_volume = xlims(2)-xlims(1);
total_num_possible_pairs = numel(C)*(numel(C)-1)/2;

rbinavg = rbinall*total_volume/(search_volume*total_num_possible_pairs); 
rbin(:,2) = rbinavg/2; 
rbin(:,3) = ubinall/pvar;

end