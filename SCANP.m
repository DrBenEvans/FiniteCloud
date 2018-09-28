function prange=SCANP(cloud,np,nb,no,ptype,rad)
%SCANP - identifies each point in the cloud's neighbouring points within a
%user-defined radius
prange=zeros(1000,np); %initialise the prange matrix
figure(3)
title('Scanning cloud for neighbours')
for ip=nb+no+1:np
  xip=cloud(1,ip);
  yip=cloud(2,ip);
  figure(3)
  scatter(cloud(1,:),cloud(2,:),'x')
  hold on
  scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
  scatter(cloud(1,nb+1:nb+no),cloud(2,nb+1:nb+no),'g')
  axis equal
  circle(xip,yip,rad)
  icount=0;
  for ich=1:np
    if(ip~=ich)
      xch=cloud(1,ich);
      ych=cloud(2,ich);
      dist=sqrt((xip-xch)^2+(yip-ych)^2);
      if(dist<rad) %check to see if line between points intersects a boundary
        [flag,inter] = CHINT([xip,yip],[xch,ych],cloud,nb,no);
        if(((ich>(nb+no)) & (flag==0)) | ((ich<=(nb+no)) & ((flag==1) & (inter==2)))) %no intersect
          icount=icount+1;
          prange(icount,ip)=ich;
          plot(xch,ych,'c+')
%        elseif((ptype==1) & (inter==1)) %a boundary point which is only intersecting with itself
%          icount=icount+1;
%          prange(icount,ip)=ich;
%          plot(xch,ych,'c+')
        end
      end
    end
   end
   hold off
 end
%
 end


