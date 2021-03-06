function [fi,dx,dy,ddx,ddy] = RBFIN2(cloud,unkno)
%RBFIN takes a point cloud and applies RBF interpolation to evaluate an
%interpolated value and compute derivatives at another point
m=2; %number of dimensions
nd=length(cloud); %number of points in cloud
xy_max=max(max(cloud));
xy_min=min(min(cloud));
r0=(xy_max-xy_min)/nd; %basis function radius
h=r0;
a1=-2.0;
a2=-1.0;
a3=0.0;
a4=1.0;
a5=2.0;
%Compute the stencil for dx, ddx derivative calculations
stencil(1,1)=cloud(1,1)+a1*h;
stencil(2,1)=cloud(2,1);
stencil(1,2)=cloud(1,1)+a2*h;
stencil(2,2)=cloud(2,1);
stencil(1,3)=cloud(1,1)+a3*h;
stencil(2,3)=cloud(2,1);
stencil(1,4)=cloud(1,1)+a4*h;
stencil(2,4)=cloud(2,1);
stencil(1,5)=cloud(1,1)+a5*h;
stencil(2,5)=cloud(2,1);
figure(5)
scatter(stencil(1,:),stencil(2,:))
hold on
scatter(cloud(1,:),cloud(2,:),'r')
axis equal
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, unkno' );
%RBF_INTERP_ND - does the RBF interpolation
fi = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 5, stencil );
%test accuracy of the interpolation
for i=1:5
    x=stencil(1,i)
    y=stencil(2,i)
    Exact=(1/sinh(pi))*sin(pi*stencil(1,i))*sinh(pi*stencil(2,i))
    Interp = fi(i)
end
%Use the stencil to compute the first and second derivatives
dx_exact=(1/sinh(pi))*pi*cos(pi*stencil(1,3))*sinh(pi*stencil(2,3))
dx=(fi(4)-fi(2))/(a4*h-a2*h)
ddx_exact=(-1/sinh(pi))*pi^2*sin(pi*stencil(1,3))*sinh(pi*stencil(2,3))
ddx_1storder=(fi(4)-2*fi(3)+fi(2))/(h*h)  %need to implement CORRECT higher order derivative approximation here...
ddx_5pointstencil=(1/(2*h*h))*(2*fi(3)-2*fi(2)+fi(1)+fi(5)-2*fi(4))
ddx_alt5pointstencil=(1/(12*h*h))*(-fi(5)+16*fi(4)-30*fi(3)+16*fi(2)-fi(1))
ddx=ddx_5pointstencil;
%Compute the stencil for dy, ddy derivative calculations
stencil(1,1)=cloud(1,1);
stencil(2,1)=cloud(2,1)+a1*h;
stencil(1,2)=cloud(1,1);
stencil(2,2)=cloud(2,1)+a2*h;
stencil(1,3)=cloud(1,1);
stencil(2,3)=cloud(2,1)+a3*h;
stencil(1,4)=cloud(1,1);
stencil(2,4)=cloud(2,1)+a4*h;
stencil(1,5)=cloud(1,1);
stencil(2,5)=cloud(2,1)+a5*h;
scatter(stencil(1,:),stencil(2,:))
hold off
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, unkno' );
%RBF_INTERP_ND - does the RBF interpolation
fi = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 5, stencil );
%test accuracy of the interpolation
for i=1:5
    x=stencil(1,i)
    y=stencil(2,i)
    Exact=(1/sinh(pi))*sin(pi*stencil(1,i))*sinh(pi*stencil(2,i))
    Interp = fi(i)
end
%Use the stencil to compute the first and second derivatives
dy_exact=(1/sinh(pi))*pi*sin(pi*stencil(1,3))*cosh(pi*stencil(2,3))
dy=(fi(4)-fi(2))/(a4*h-a2*h)
ddy_exact=(1/sinh(pi))*pi^2*sin(pi*stencil(1,3))*sinh(pi*stencil(2,3))
ddy_1storder=(fi(4)-2*fi(3)+fi(2))/(h*h)  %need to implement CORRECT higher order derivative approximation here...
ddy_5pointstencil=(1/(2*h*h))*(2*fi(3)-2*fi(2)+fi(1)+fi(5)-2*fi(4))
ddy_alt5pointstencil=(1/(12*h*h))*(-fi(5)+16*fi(4)-30*fi(3)+16*fi(2)-fi(1))
ddy=ddy_5pointstencil;
end