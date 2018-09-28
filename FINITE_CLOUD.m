%% FINITE CLOUD TEST CODE
%% this MATLAB program is designed to solve Laplace's equation around an arbitrary shape
%% defined by a set of points (in a .prf file) using the 'Finite Cloud' numerical scheme
%%
%% written by B. Evans (Swansea University) - coding started on 27th May 2016... continuing summer 2018...

%% First clear variables/figures and load in the geometry:
clc;clear;close all;
%% User input variables declared here:
filename = 'circle.prf'
ptins = 0.001;  %threshold for Levy flight point insertion
nLevy=5; %number of Levy flight steps
beta=1.5;  %Levy flight Beta parameter (set betweeen 1.0 and 2.0)
bstep=0.1;  %step size away from boundary for points inserted from the boundary
tolp=0.001; %the smallest distance allowed between points in the cloud
rad=1.0;  %the 'influence radius' around each point
gradlim=1;  %gradient limiter
dt=0.1;  %timestep
nstep=100; %number of timesteps
mxpnt=1000; %maximum allowed number of points
%% First read in the geometry and create the starting point master point cloud
[cloud,np,nb,no,ptype,pinsert] = RGEOM(filename); 
%% Then check for duplicate points
[cloud,np,nb,no,ptype,pinsert] = CHDUP(cloud,np,nb,no,ptype,pinsert,tolp);
%% next compute the normals of all the points on the boundary
[pnorm] = CMPNO(np,nb,no,cloud);
%% display the initial point cloud
figure(1)
scatter(cloud(1,:),cloud(2,:),'x')
axis equal
title('Initial Point Cloud Distribution and Boundary Normals')
hold on
quiver(cloud(1,1:(nb+no)),cloud(2,1:(nb+no)),pnorm(1,:),pnorm(2,:))
%% perform initial point insertion
[cloud,np,ptype,pinsert] = PNTIN(cloud,np,nb,no,ptype,pinsert,pnorm,nLevy,beta,bstep);
%% TEMPORARY OVERWRITE OF CLOUD AND UNKNO FOR CODE TESTING %%%%%%%%%%%%%%%
clear cloud
clear unkno
clear ptype
coord=linspace(0,1,20);
ip=0;
for ix=1:20
for iy=1:20
ip=ip+1;
cloud(ip,1)=coord(ix);
cloud(ip,2)=coord(iy);
end
end
%add in a surrounding 'buffer zone'
np=length(cloud);
ip=np;
for i=1:20
    ip=ip+1;
    cloud(ip,1)=coord(i);
    cloud(ip,2)=-0.05;
end
for i=1:20
    ip=ip+1;
    cloud(ip,1)=coord(i);
    cloud(ip,2)=1.05;
end
for i=1:20
    ip=ip+1;
    cloud(ip,1)=-0.05;
    cloud(ip,2)=coord(i);
end
for i=1:20
    ip=ip+1;
    cloud(ip,1)=1.05;
    cloud(ip,2)=coord(i);
end
cloud=transpose(cloud);
scatter(cloud(1,:),cloud(2,:))
np=length(cloud)
for ip=1:np
    %if((cloud(2,ip)>0.1)&(cloud(2,ip)<0.9)&(cloud(1,ip)>0.1)&(cloud(1,ip)<0.9))
    %    unkno(ip)=(1/sinh(pi))*sin(pi*cloud(1,ip))*sinh(pi*cloud(2,ip))+sin(pi*cloud(2,ip));
    %else
        unkno(ip)=(1/sinh(pi))*sin(pi*cloud(1,ip))*sinh(pi*cloud(2,ip));
    %end
end
ptype(1:np)=3;
for ip=1:(np-80)
    if((cloud(1,ip)==0)|(cloud(2,ip)==0)|(cloud(1,ip)==1)|(cloud(2,ip)==1)) %domain edge
       ptype(ip)=2;
    end
end
for ip=401:np
    ptype(ip)=2;
end
nb=0;
no=0;
%% display the new point cloud
figure(2)
scatter(cloud(1,:),cloud(2,:),'x')
hold on
scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
scatter(cloud(1,nb+1:nb+no),cloud(2,nb+1:nb+no),'g')
axis equal
title('Iteration 1 Point Cloud Distribution')
%% scan the point cloud to identify which other points are in range of influence
prange=SCANP(cloud,np,nb,no,ptype,rad);
%% Apply initial condition to cloud points
%unkno=INCON(cloud,ptype,np,nb);
%% display updated solution
   figure(7)
   scatter(cloud(1,:),cloud(2,:),'x')
   hold on
   scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
   axis equal
   title('Updated Point Cloud Distribution and unknown distribution')
   % visualise the unknown
   tri=delaunay(cloud(1,:),cloud(2,:));
   figure(7)
   trisurf(tri,cloud(1,:),cloud(2,:),unkno,unkno)
   colorbar
   %caxis([-4,4])
   view(0,90) 
%%  TIMESTEP, UPDATE NODAL VALUES AND APPLY BOUNDARY CONDITION
for it=1:nstep
    fprintf('************************\n');
    fprintf('Iteration number: %d\n',it);
    fprintf('************************\n\n');
    npstore(it)=np;
   %% USE RADIAL BASIS FUNCTION INTERP (RBFIN) TO COMPUTE THE LAPLACIAN OPERATOR AT EACH POINT
   [residual,resnorm]=RESID(cloud,unkno,ptype,prange,np,nb,gradlim);
   figure(8)
   istep(it)=it;
   resplot(it)=resnorm;
   plot(istep,resplot);
   xlabel('Timestep')
   ylabel('Residual Norm')
   % visualise the residual distribution
   tri=delaunay(cloud(1,:),cloud(2,:));
   figure(6)
   trisurf(tri,cloud(1,:),cloud(2,:),residual,residual)
   title('Delaunay/trisurf representation of the residual')
   view(0,90)
   hold on
   scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
   axis equal
   colorbar
   % update the solution (timestep)
   unkno=TIMES(residual,np,nb,unkno,ptype,dt); 
   %% display updated solution
   figure(7)
   scatter(cloud(1,:),cloud(2,:),'x')
   hold on
   scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
   axis equal
   title('Updated Point Cloud Distribution and unknown distribution')
   % visualise the unknown
   tri=delaunay(cloud(1,:),cloud(2,:));
   figure(7)
   trisurf(tri,cloud(1,:),cloud(2,:),unkno,unkno)
   colorbar
   %caxis([-4,4])
   view(0,90) 
   %% Test for point insertion and update pinsert
%   [pinsert]=INSERT(np,nb,residual,pinsert,ptins);
   %% perform point insertion if np<mxpnt
%   if(np<mxpnt)
%     [cloud,np,ptype,pinsert] = PNTIN(cloud,np,nb,no,ptype,pinsert,pnorm,nLevy,beta,bstep);
     %% Then check for duplicate points
%     [cloud,np,nb,no,ptype,pinsert] = CHDUP(cloud,np,nb,no,ptype,pinsert,tolp);
     %% scan the point cloud to identify which other points are in range of influence
%     prange=SCANP(cloud,np,nb,no,ptype,rad);
     %% Interpolate unknowns to the newly inserted points
%     unkno = UNKUP(unkno,cloud,np,prange);
%   end
   %% display the new point cloud
%   figure(7)
%   scatter(cloud(1,:),cloud(2,:),'x')
%   hold on
%   scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
   %axis equal
%   title('Updated Point Cloud Distribution and unknown distribution')
   % visualise the unknown
%   tri=delaunay(cloud(1,:),cloud(2,:));
%   figure(7)
%   trisurf(tri,cloud(1,:),cloud(2,:),unkno,unkno)
%   colorbar
%   caxis([-4,4])
   %view(0,90) 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
