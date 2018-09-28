function [unkno] = UNKUP(unkno,cloud,np,prange)
%UNKUP updates the unkno vector (of unknowns) after new cloud point insertion
orignp=length(unkno);
for ip=orignp+1:np
    icount=1;
    unkno_loc(1)=0.0;
    cloud_loc(:,1)=cloud(:,ip);
    for icheck=1:1000
        ip2=prange(icheck,ip);
        if((ip2>0) & (ip2<=orignp))
          icount=icount+1;
          cloud_loc(:,icount)=cloud(:,ip2);
          unkno_loc(:,icount)=unkno(ip2);
        end
    end
    
    %check to see if there were points within range
    if(length(unkno_loc)<2)  %this implies that there are no points within range
       unkno(ip)=0.0;
    else
       unkno_loc(1)=mean(unkno_loc(2:icount));  %require an initial 'guess' at the value at the new point
       %call RBFIN to do the interpolation and compute second derivatives
      [fi,dx,dy,ddx,ddy] = RBFIN(cloud_loc,unkno_loc);
      unkno(ip)=fi(9);
    end
end

