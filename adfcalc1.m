function [rbin,rvec,theta]=adfcalc1(C,rvec,theta,xlims,ylims,varargin)
% Parameters:
%   C       matrix (n,2) OR (n,3) of n (x,y) OR (x,y,z) particle locations
%   rvec    row vector of radius segments for particles
%   theta   row vector of angle segments for particles (in radians) (theta =
%       angle from vertical axis)
%   xlim    vector of min,max x location
%   ylim    vector of min,max y location
%   zlim    vector of min,max z location (OPTIONAL)
%
% This function takes as input a matrix containing the co-ordinates of the
% particles, a vector of shell radii, and a vector of angles, and finds the
% ADF. 
%
% Note: theta=0 means P2 is directly above or below P1. 
%       theta=pi/2 means P2 is in a horizontal plane with P1.
%       theta is defined for [0,pi/2] radians.
%
% rbin(:,i) is the ADF vector, normalized by the ADF of a uniform
% distribution, for the i-th theta bin.

if size(C,2) == 2
    dim2 = true;
else
    dim2 = false;
    zlims = varargin{1};
end

k = size(theta,2);
m = size(rvec,2);

rbin = zeros(m-1,k-1);

if dim2
    % subselect set of first particles (P1)
    subselect1 = (C(:,1)>(xlims(1)+rvec(end)) & C(:,1)<(xlims(2)-rvec(end)) & C(:,2)<(ylims(2)-rvec(end)) & C(:,2)>(ylims(1)+rvec(end)));
    % subselect set of second particles (P2)
    subselect2 = (C(:,1)>(xlims(1)) & C(:,1)<(xlims(2)) & C(:,2)<(ylims(2)) & C(:,2)>(ylims(1)));
else
    subselect1 = (C(:,1)>(xlims(1)+rvec(end)) & C(:,1)<(xlims(2)-rvec(end)) & C(:,2)<(ylims(2)-rvec(end)) & C(:,2)>(ylims(1)+rvec(end)) & C(:,3)<(zlims(2)-rvec(end)) & C(:,3)>(zlims(1)+rvec(end)));
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
clear subselect1
clear subselect2
n1 = length(x1)
n2 = length(x2)
% r12 = zeros(n1,n2);
% th12 = zeros(n1,n2);

% precompute distances and angles
if dim2
    for i = 1:n1   % for all P1
        % find distance to P2
        r12 = sqrt((x1(i)-x2(:)).^2 + (y1(i)-y2(:)).^2);
        % find angle theta to P2
        th12 = abs(atan((x2(:)-x1(i))./(y2(:)-y1(i))));
        for j = 1:m-1   % for each annulus
            for h = 1:k-1   % for each angle bin
                rbin(j,h) = rbin(j,h) + sum(r12 >= rvec(j) & r12 <= rvec(j+1) & th12 >= theta(h) & th12 <= theta(h+1));
            end
        end
    end
else
    for i = 1:n1   % for all P1
%         if ~mod(i,round(n1/5))
%             fprintf(['i = ' num2str(i) ' of ' num2str(n1) '\n']);
%         end
        % find distance to P2
        r12 = sqrt((x1(i)-x2(:)).^2 + (y1(i)-y2(:)).^2 + (z1(i)-z2(:)).^2);
        % find angle theta to P2
        th12 = acos(abs(z1(i)-z2(:))./r12);
        for j = 1:m-1   % for each annulus
            for h = 1:k-1   % for each angle
                rbin(j,h) = rbin(j,h) + sum(r12 >= rvec(j) & r12 < rvec(j+1) & th12 >= theta(h) & th12 < theta(h+1));
            end
        end
    end
end

% normalization values
if dim2
	avg_pd = n2/(diff(xlims)*diff(ylims));  % average particle density
else
	avg_pd = n2/(diff(xlims)*diff(ylims)*diff(zlims));  % average particle density
end
    
area1=zeros(m-1,k-1);
for j = 1:m-1     % for each annulus
    for i = 1:k-1   % for each angle
        % area or volume of bin
        if dim2
            area1(j,i) = pi*((rvec(j+1))^2-rvec(j)^2)*(theta(i+1)-theta(i))/(2*pi);
        else
            area1(j,i) = ( ( (Vcone(theta(i+1),rvec(j+1)) + Vcap(theta(i+1),rvec(j+1))) - ...
                (Vcone(theta(i),rvec(j+1)) + Vcap(theta(i),rvec(j+1))) ) - ...
               ( (Vcone(theta(i+1),rvec(j)) + Vcap(theta(i+1),rvec(j))) - ...
                (Vcone(theta(i),rvec(j)) + Vcap(theta(i),rvec(j))) ) )*2;
        end
        % find radius and angle bin of P2s
        % normalize by number expected if distribution were random, averaged over P1 locations
        rbin(j,i) = rbin(j,i)/(area1(j,i)*avg_pd*n1);
    end
end

rvec = rvec(2:end)-0.5*diff(rvec); % radius bin centers
theta = theta(2:end)-0.5*diff(theta); % angle bin centers

end

function V = Vcone(th,r)
% volume of a cone
V = pi/3*r^3*sin(th)^2*cos(th);
end

function V = Vcap(th,r)
% volume of a spherical cap
V = pi/3*r^3*(2 - 3*cos(th) + cos(th)^3);
end