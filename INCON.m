function [unkno]=INCON(cloud,ptype,np,nb)
%INCON - sets up the initial condition unknowns on all cloud points
special=1;
if(special==1)  %special case solving the w=w_o*sin(pi*x) upper boundary scenario -> with known analytical solution
    coord=linspace(0,1,20);
    w_0=1.0;
    for i=1:18
        unkno(i)=0.0;
    end
    for i=19:38
        unkno(i)=w_0*sin(pi*coord(i-18));
    end
    for i=39:76
        unkno(i)=0.0;
    end
    for i=77:np
        unkno(i)=0.0;  %initial condition in the domain
    end 
else
  for i=1:np
      if(ptype(i)==1)  %shape boundary
          unkno(i)=0;
      elseif(ptype(i)==2) %outer domain
%          unkno(i)=0;
           x=cloud(1,i);
           y=cloud(2,i);
           if(x>=0)
             theta=atan(y/x);
           else
             theta=pi+atan(y/x);
           end
           unkno(i)=4*sin(5*theta);
      elseif(ptype(i)==3)  %domain points
          unkno(i)=0;
      end
  end
end
%visualise this inital condition set up via Delaunay triangulation
tri=delaunay(cloud(1,:),cloud(2,:));
figure(4)
trisurf(tri,cloud(1,:),cloud(2,:),unkno,unkno)
title('Delaunay/trisurf representation of the initial condition')
view(0,90)
hold on
scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
axis equal
colorbar
end

