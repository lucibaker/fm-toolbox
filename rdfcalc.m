function [rbin]=rdfcalc(C,rvecarray,xlims,ylims,varargin)
%
% Parameters:
%   C       matrix (n,2) OR (n,3) of n (x,y) OR (x,y,z) particle locations
%   rvec    array of radi segments for particles
%   xlim    vector of min,max x location
%   ylim    vector of min,max y location
%   zlim    vector of min,max z location (OPTIONAL)
%
% This function takes as input a matrix containing the co-ordinates of the
% particles and a vector of shell radii, and finds the RDF
%
% rbin(:,1) is the radius vector in meters. rbin(:,2) is the RDF
% vector, normalized by the RDF of a uniform distribution.

if size(C,2) == 2
    dim2 = true;
else
    dim2 = false;
    zlims = varargin{1};
end

m = size(rvecarray,2);
rbin = zeros(m-1,2);
rbin(:,1)=rvecarray(2:end)-0.5*diff(rvecarray);
if dim2
    % subselect set of first particles (P1)
    subselect1 = (C(:,1)>(xlims(1)+rvecarray(end)) & C(:,1)<(xlims(2)-rvecarray(end)) & C(:,2)<(ylims(2)-rvecarray(end)) & C(:,2)>(ylims(1)+rvecarray(end)));
    % subselect set of second particles (P2)
    subselect2 = (C(:,1)>(xlims(1)) & C(:,1)<(xlims(2)) & C(:,2)<(ylims(2)) & C(:,2)>(ylims(1)));
else
    subselect1 = (C(:,1)>(xlims(1)+rvecarray(end)) & C(:,1)<(xlims(2)-rvecarray(end)) & C(:,2)<(ylims(2)-rvecarray(end)) & C(:,2)>(ylims(1)+rvecarray(end)) & C(:,3)<(zlims(2)-rvecarray(end)) & C(:,3)>(zlims(1)+rvecarray(end)));
    subselect2 = (C(:,1)>(xlims(1)) & C(:,1)<(xlims(2)) & C(:,2)<(ylims(2)) & C(:,2)>(ylims(1)) & C(:,3)<(zlims(2)) & C(:,3)>(zlims(1)));
end
x1 = C(subselect1,1);
x2 = C(subselect2,1);
y1 = C(subselect1,2);
y2 = C(subselect2,2);
if ~dim2
    z1 = C(subselect1,3);
    z2 = C(subselect2,3);
end
clear C
n1 = length(x1);
n2 = length(x2);
% r12 = zeros(1,n2);

% precompute distances
if dim2
    for i = 1:n1   % for all P1
        % find distance to P2
        r12 = sqrt((x1(i)-x2(:)).^2 + (y1(i)-y2(:)).^2);
        for j = 1:m-1   % count for each annulus
            rbin(j,2) = rbin(j,2) + sum(sum(r12 >= rvecarray(j) & r12 <= rvecarray(j+1)));
        end
    end
else
    for i = 1:n1   % for all P1
%         if ~mod(i,round(n1/10))
%             fprintf(['i = ' num2str(i) ' of ' num2str(n1) '\n']);
%         end
        % find distance to P2
        r12 = sqrt((x1(i)-x2(:)).^2 + (y1(i)-y2(:)).^2 + (z1(i)-z2(:)).^2);
        for j = 1:m-1   % count for each annulus
            rbin(j,2) = rbin(j,2) + sum(sum(r12 >= rvecarray(j) & r12 <= rvecarray(j+1)));
        end
    end
end

% normalization values
if dim2
	avg_pd = n2/(diff(xlims)*diff(ylims));  % average particle density
else
	avg_pd = n2/(diff(xlims)*diff(ylims)*diff(zlims));  % average particle density
end


for j = 1:m-1     % for each annulus
    % area or volume of bin
    if dim2
        area1 = pi*((rvecarray(j+1))^2-rvecarray(j)^2);
    else
        area1 = 4/3*pi*((rvecarray(j+1))^3-rvecarray(j)^3);
    end
    % find radius bin of P2s
    % normalize by number expected if distribution were random, averaged over P1 locations
    rbin(j,2) = rbin(j,2)/(area1*avg_pd*n1);
end

end
