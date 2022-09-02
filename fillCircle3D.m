function points = fillCircle3D(center,normal,radius,c)
% center, normal can be column vectors
N = size(center,1);
theta = linspace(0,2*pi,20);
for i = 1:N
    v=null(normal(i,:));
    points=repmat(center(i,:)',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    if nargout == 0
        fill3(points(1,:),points(2,:),points(3,:),c); hold on
    end
end
end