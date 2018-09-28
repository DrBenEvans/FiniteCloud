function [residual,resnorm]=RESID(cloud,unkno,ptype,prange,np,nb,gradlim)
%RESIDUAL - computes the point residuals for the governing equation (in
%this case Laplace) by radial basis function interpolation
residual=zeros(1,np);
for ip=1:np
    clear cloud_loc
    clear unkno_loc
    if(ptype(ip)==3)
      %construct the local cloud of points (within range)
      icount=1;
      cloud_loc(:,1)=cloud(:,ip);
      unkno_loc(1)=unkno(ip);
      for icheck=1:1000
          ip2=prange(icheck,ip);
          if(ip2>0)
              icount=icount+1;
              cloud_loc(:,icount)=cloud(:,ip2);
              unkno_loc(:,icount)=unkno(ip2);
          end
      end
      %check to see if there were points within range
      if(length(unkno_loc)<2)  %this implies that there are no points within range
          residual(ip)=0;
      else
         %call rbfDerivative to do the interpolation and compute second derivatives
         [fi,dx,dy,ddx,ddy]=rbfDerivative(cloud_loc,unkno_loc);
         %calculate the residual at each point (applying the limiter)
         if((abs(dx)>gradlim)|(abs(dy)>gradlim)|(abs(ddx)>gradlim)|(abs(ddy)>gradlim))
             residual(ip)=gradlim;
          if((ddx+ddy)>gradlim)
%              residual(ip)=gradlim;
               residual(ip)=0;
          elseif((ddx+ddy)<-gradlim)
%              residual(ip)=-gradlim;
               residual(ip)=0;
          elseif((isnan(dx)==1)|(isnan(dy)==1)|(isnan(ddx)==1)|(isnan(ddy)==1))
             residual(ip)=0;
          else
             residual(ip)=ddx+ddy;  %LAPLACE EQUATION!!
          end
      end
    end
end
%compute the mean of the residual as a measure of convergence
resnorm=mean(abs(residual));
end

