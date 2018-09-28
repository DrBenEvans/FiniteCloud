clear all
close all
clc
coord=linspace(0,1,20);
ip=0;
for ix=1:20
for iy=1:20
ip=ip+1;
cloud(ip,1)=coord(ix);
cloud(ip,2)=coord(iy);
end
end
cloud=transpose(cloud);
scatter(cloud(1,:),cloud(2,:))
np=length(cloud)
for ip=1:np
unkno(ip)=(1/sinh(pi))*sin(pi*cloud(1,ip))*sinh(pi*cloud(2,ip));
end
% visualise the residual distribution
  tri=delaunay(cloud(1,:),cloud(2,:));
  figure
  trisurf(tri,cloud(1,:),cloud(2,:),unkno,unkno)
  title('Delaunay/trisurf representation of the unknown distribution')
  view(0,90)
  axis equal
  colorbar