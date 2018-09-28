function [cloud,np,nb,no,ptype,pinsert] = RGEOM(filename)
%RGEOM reads the shape boundary points and sets up initial arrays cloud, ptype and pinsert
%N.B. ptype=1 - internal shape boundary, ptype=2 - outer boundary, ptype=3 - internal point
%  
initialgeom=3; %flag to dictate the type of outer point distribution
%
fid = fopen(filename);
cloud = fscanf(fid,'%f',[2,inf]);
fclose(fid);
np=length(cloud);
nb=np; % number of shape boundary points
%create the point identifier vector and point insertion check vector
for ip=1:np
    ptype(ip)=1;  %shape boundary point
    pinsert(ip)=1;   %1=insert new point
end
%identify rough size of geometry
max_xy=max(max(cloud));
min_xy=min(min(cloud));
scale=max_xy-min_xy;
%identify the rough centre of the shape
av_x=mean(cloud(1,:));
av_y=mean(cloud(2,:));
if(initialgeom==1)
%add outer boundary points as four corner points
  cloud(1,np+1)=av_x+scale*10;
  cloud(2,np+1)=av_y+scale*10;
  ptype(np+1)=2; %outer boundary
  pinsert(np+1)=0;  %do not insert
  cloud(1,np+2)=av_x+scale*10;
  cloud(2,np+2)=av_y-scale*10;
  ptype(np+2)=2; %outer boundary
  pinsert(np+2)=0;  %do not insert
  cloud(1,np+3)=av_x-scale*10;
  cloud(2,np+3)=av_y-scale*10;
  ptype(np+3)=2; %outer boundary
  pinsert(np+3)=0;  %do not insert
  cloud(1,np+4)=av_x-scale*10;
  cloud(2,np+4)=av_y+scale*10;
  ptype(np+4)=2; %outer boundary
  pinsert(np+4)=0;  %do not insert
  %new number of points in the cloud
  np=length(cloud);
  no=4; %number of outer boundary points
elseif(initialgeom==2)
    %add a circular outer boundary of radius 10.0 (laplace solution test
    %case)
    %set value of geometric parameters
    a=0.0;  %x-coordinate of circle centre
    b=0.0;  %y-coordiante of circle centre
    c=10.0;  %circle radius
    %create the circle points (on  2 curves of 50 points each)
    c1=linspace(0,-pi,50);
    c2=linspace(0,pi,50);
    x1=c*sin(c1);
    y1=-c*cos(c1);
    x2=c*sin(c2);
    y2=c*cos(c2);
    x1=x1+a;
    x2=x2+a;
    y1=y1+b;
    y2=y2+b;     
    for i=1:50
      cloud(1,np+i)=x1(i);
      cloud(2,np+i)=y1(i);
      pinsert(np+i)=1; %insert points from here
      ptype(np+i)=2; %outer boundary
    end
    for i=2:49
      cloud(1,np+49+i)=x2(i);
      cloud(2,np+49+i)=y2(i);
      pinsert(np+49+i)=1; %insert points from here
      ptype(np+49+i)=2; %outer boundary
    end
    %new number of points in the cloud
    np=length(cloud);
    no=98; %number of outer boundary points
elseif(initialgeom==3)
    cloud=0.0;
    %rectangular domain with no inner body
    coord=linspace(0,1,20);
    ip=0;
    for ix=1:20
        ip=ip+1;
        cloud(ip,1)=coord(ix);
        cloud(ip,2)=0.0;
        pinsert(ip)=1;
        ptype(ip)=2;
    end
    for iy=2:19
        ip=ip+1;
        cloud(ip,1)=1.0;
        cloud(ip,2)=coord(iy);
        pinsert(ip)=1;
        ptype(ip)=2;
    end
    for ix=1:20
        ip=ip+1;
        cloud(ip,1)=coord(21-ix);
        cloud(ip,2)=1.0;
        pinsert(ip)=1;
        ptype(ip)=2;
    end
    for iy=2:19
        ip=ip+1;
        cloud(ip,1)=0.0;
        cloud(ip,2)=coord(21-iy);
        pinsert(ip)=1;
        ptype(ip)=2;
    end
    cloud=transpose(cloud);
    no=ip;
    nb=0;
    np=length(cloud);
    %reverse order
    for ip=1:np
        cloud_new(1,ip)=cloud(1,np+1-ip);
        cloud_new(2,ip)=cloud(2,np+1-ip);
    end
    cloud=cloud_new;
end


