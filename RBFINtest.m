function [fi,dx,dy,ddx,ddy] = RBFINtest(point,cloud,unkno)
%RBFIN takes a point cloud and applies RBF interpolation to evaluate an
%interpolated value and compute derivatives at another point
%TESTING AND DEBUGGING VERSION ONLY!!
m=2; %number of dimensions
nd=length(cloud); %number of points in cloud
xy_max=max(max(cloud));
xy_min=min(min(cloud));
r0=(xy_max-xy_min)/nd; %basis function radius
stensize=20*r0;
%Compute the stencil for derivative calculations
stencil(1,1)=point(1)-2*stensize;
stencil(1,2)=point(2);
stencil(2,1)=point(1)-stensize;
stencil(2,2)=point(2);
stencil(3,1)=point(1)+stensize;
stencil(3,2)=point(2);
stencil(4,1)=point(1)+2*stensize;
stencil(4,2)=point(2);
stencil(5,1)=point(1);
stencil(5,2)=point(2)-2*stensize;
stencil(6,1)=point(1);
stencil(6,2)=point(2)-stensize;
stencil(7,1)=point(1);
stencil(7,2)=point(2)+stensize;
stencil(8,1)=point(1);
stencil(8,2)=point(2)+2*stensize;
stencil(9,1)=point(1);
stencil(9,2)=point(2);
stencil=stencil';
figure(5)
scatter(stencil(1,:),stencil(2,:))
hold on
scatter(cloud(1,:),cloud(2,:),'r')
hold off
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, unkno' );
%RBF_INTERP_ND - does the RBF interpolation
fi = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 9, stencil );
%Use the stencil to compute the first and second derivatives
dx=(fi(3)-fi(2))/(2*stensize);
dy=(fi(7)-fi(6))/(2*stensize);
dxright=(fi(4)-fi(9))/(2*stensize);
dxleft=(fi(9)-fi(1))/(2*stensize);
ddx=(dxright-dxleft)/(2*stensize);
%ddx=(fi(4)-2*fi(9)-fi(1))/(2*r0)^2;
dytop=(fi(8)-fi(9))/(2*stensize);
dybot=(fi(9)-fi(5))/(2*stensize);
ddy=(dytop-dybot)/(2*stensize);
%ddy=(fi(8)-2*fi(9)-fi(5))/(2*r0)^2;
%figure(5)
%tri = delaunay ( cloud(1,:), cloud(2,:) );
%trisurf ( tri, cloud(1,:), cloud(2,:), unkno );
%  xlabel ( '<--- X --->' );
%  ylabel ( '<--- Y --->' );
%  zlabel ( '<---Z(X,Y)--->' );
%  title ( 'Exact function' )
end