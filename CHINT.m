function [flag,inter] = CHINT(L1,L2,cloud,nb,no)
%CHINT - this routine checks whether a Levy flight step (defined by points
%L1 and L2) intersects with any of the shape boundary (defined by points B1 and B2) and returns 1 for yes
%(intersect) and 0 for no (no intersect)
% ALGORITHM by D. Naumann
flag=0;
inter=0;
%check for intersection with shape boundary
for ip=1:nb
  B1(1)=cloud(1,ip);
  B1(2)=cloud(2,ip);
  if(ip==nb)
   B2(1)=cloud(1,1);
   B2(2)=cloud(2,1);
  else
   B2(1)=cloud(1,ip+1);
   B2(2)=cloud(2,ip+1);
  end
  
  % calculate the intersection point
  A = B2(2)-B1(2);
  B = B1(1)-B2(1);
  C = A*B1(1) + B*B1(2);

  Atarget = L2(2)-L1(2);
  Btarget = L1(1)-L2(1);
  Ctarget = Atarget*L1(1) + Btarget*L1(2);

  detAB = A*Btarget - Atarget*B;
  detCB = C*Btarget - Ctarget*B;
  detAC = A*Ctarget - Atarget*C;

  if(abs(detAB) < abs(10e-12))  %lines are parallel
     xtemp=-1000;
     ytemp=-1000; 
  else %non-parallel
     xtemp=detCB/detAB;
     ytemp=detAC/detAB;
  end

  %check if intersection lies on the boundary segment considered
  if(((xtemp - min([B1(1),B2(1)]))>-10e-12) & ((xtemp - max([B1(1),B2(1)]))<10e-12) & ((ytemp - min([B1(2),B2(2)]))>-10e-12) & ((ytemp - max([B1(2),B2(2)]))<10e-12)) 
    if(((xtemp - min([L1(1),L2(1)]))>-10e-12) & ((xtemp - max([L1(1),L2(1)]))<10e-12) & ((ytemp - min([L1(2),L2(2)]))>-10e-12) & ((ytemp - max([L1(2),L2(2)]))<10e-12))
       intersect(1)=xtemp;
       intersect(2)=ytemp;
       flag=1;
       inter=inter+1;
    end
  end
end
%
%check for intersection with outer boundary
for ip=nb+1:nb+no
  B1(1)=cloud(1,ip);
  B1(2)=cloud(2,ip);
  if(ip==(nb+no))
   B2(1)=cloud(1,nb+1);
   B2(2)=cloud(2,nb+1);
  else
   B2(1)=cloud(1,ip+1);
   B2(2)=cloud(2,ip+1);
  end
  
  % calculate the intersection point
  A = B2(2)-B1(2);
  B = B1(1)-B2(1);
  C = A*B1(1) + B*B1(2);

  Atarget = L2(2)-L1(2);
  Btarget = L1(1)-L2(1);
  Ctarget = Atarget*L1(1) + Btarget*L1(2);

  detAB = A*Btarget - Atarget*B;
  detCB = C*Btarget - Ctarget*B;
  detAC = A*Ctarget - Atarget*C;

  if(abs(detAB) < abs(10e-12))  %lines are parallel
     xtemp=-1000;
     ytemp=-1000; 
  else %non-parallel
     xtemp=detCB/detAB;
     ytemp=detAC/detAB;
  end

  %check if intersection lies on the boundary segment considered
  if(((xtemp - min([B1(1),B2(1)]))>-10e-12) & ((xtemp - max([B1(1),B2(1)]))<10e-12) & ((ytemp - min([B1(2),B2(2)]))>-10e-12) & ((ytemp - max([B1(2),B2(2)]))<10e-12)) 
    if(((xtemp - min([L1(1),L2(1)]))>-10e-12) & ((xtemp - max([L1(1),L2(1)]))<10e-12) & ((ytemp - min([L1(2),L2(2)]))>-10e-12) & ((ytemp - max([L1(2),L2(2)]))<10e-12))
       intersect(1)=xtemp;
       intersect(2)=ytemp;
       flag=1;
       inter=inter+1;
    end
  end
end

end

