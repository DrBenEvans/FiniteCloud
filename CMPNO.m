function [pnorm] = CMPNO(np,nb,no,cloud);
%CMPNO - computes the outward-facing unit normals of all points on the
%boundary
if(nb>0)
  pnorm=zeros(2,nb);
  %first point on shape boundary
  xleft=cloud(1,nb);
  yleft=cloud(2,nb);
  xright=cloud(1,2);
  yright=cloud(2,2);
  tang(1)=xleft-xright;
  tang(2)=yleft-yright;
  pnorm(1,1)=-tang(2);
  pnorm(2,1)=tang(1);
  for ip=2:nb-1
    xleft=cloud(1,ip-1);
    yleft=cloud(2,ip-1);
    xright=cloud(1,ip+1);
    yright=cloud(2,ip+1);
    tang(1)=xleft-xright;
    tang(2)=yleft-yright;
    pnorm(1,ip)=-tang(2);
    pnorm(2,ip)=tang(1);
  end
  %last point on boundary
  xleft=cloud(1,nb-1);
  yleft=cloud(2,nb-1);
  xright=cloud(1,1);
  yright=cloud(2,1);
  tang(1)=xleft-xright;
  tang(2)=yleft-yright;
  pnorm(1,nb)=-tang(2);
  pnorm(2,nb)=tang(1);
  %normalise
  for ip=1:nb
    magno=sqrt((pnorm(1,ip)^2+pnorm(2,ip)^2));
    pnorm(1,ip)=pnorm(1,ip)/magno;
    pnorm(2,ip)=pnorm(2,ip)/magno;
  end
end
%first point on the outer boundary
xleft=cloud(1,nb+no);
yleft=cloud(2,nb+no);
xright=cloud(1,nb+2);
yright=cloud(2,nb+2);
tang(1)=xleft-xright;
tang(2)=yleft-yright;
pnorm(1,1+nb)=-tang(2);
pnorm(2,1+nb)=tang(1);
for ip=nb+2:nb+no-1
    xleft=cloud(1,ip-1);
    yleft=cloud(2,ip-1);
    xright=cloud(1,ip+1);
    yright=cloud(2,ip+1);
    tang(1)=xleft-xright;
    tang(2)=yleft-yright;
    pnorm(1,ip)=-tang(2);
    pnorm(2,ip)=tang(1);
end
%last point on outer boundary
xleft=cloud(1,nb+no-1);
yleft=cloud(2,nb+no-1);
xright=cloud(1,nb+1);
yright=cloud(2,nb+1);
tang(1)=xleft-xright;
tang(2)=yleft-yright;
pnorm(1,nb+no)=-tang(2);
pnorm(2,nb+no)=tang(1);
%normalise
for ip=nb+1:nb+no
    magno=sqrt((pnorm(1,ip)^2+pnorm(2,ip)^2));
    pnorm(1,ip)=pnorm(1,ip)/magno;
    pnorm(2,ip)=pnorm(2,ip)/magno;
end

