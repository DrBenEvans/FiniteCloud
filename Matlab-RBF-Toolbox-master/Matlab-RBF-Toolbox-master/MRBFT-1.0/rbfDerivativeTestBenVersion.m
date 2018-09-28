function [fi,dx,dy,ddx,ddy]=rbfDerivativeTestBenVersion(point,cloud,unkno)
% this function takes in a point cloud and unknown matrix, generates the
% derivative distribution across the cloud and interpolates dx, dy, ddx,
% ddy and a given point
% ----------------------------------------------------------------------
% Input parameters and initiate RBF
phi = gax();
s = 4.5; 
nCh = inf;                  % norm choice
G = F2c;

%[x, y] = rbfCenters.circleCenters(N,true,1,1,false);
%f = G.F(x,y);

x=cloud(1,:);
x=transpose(x);
y=cloud(2,:);
y=transpose(y);
f=transpose(unkno);

[r, rx, ry] = phi.distanceMatrix2d(x,y);
B = phi.rbf(r,s);

mu = 0;
safe = true;  % use mldivide rather than Cholesky directly
a = rbfx.solve(B,f,mu,safe);


H = phi.D1(r,s,rx);
fx = H*a;
%fprintf('dx error = %4.2e\n',norm(fx - G.x1(x,y), nCh));

H = phi.D2(r,s,rx);
fxx = H*a;
%fprintf('dxx error = %4.2e\n',norm(fxx - G.x2(x,y), nCh));

H = phi.D3(r,s,rx);
fxxx = H*a;
%fprintf('dxxx error = %4.2e\n',norm(fxxx - G.x3(x,y), nCh));

H = phi.D4(r,s,rx);
fxxxx = H*a;
%fprintf('dxxxx error = %4.2e\n',norm(fxxxx - G.x4(x,y), nCh));

H = phi.D1(r,s,ry);
fy = H*a;
%fprintf('dy error = %4.2e\n',norm(fy - G.y1(x,y), nCh));

H = phi.D2(r,s,ry);
fyy = H*a;
%fprintf('dyy error = %4.2e\n',norm(fyy - G.y2(x,y), nCh));

H = phi.D3(r,s,ry);
fyyy = H*a;
%fprintf('dyyy error = %4.2e\n',norm(fyyy - G.y3(x,y), nCh));

H = phi.D4(r,s,ry);
fyyyy = H*a;
%fprintf('dyyyy error = %4.2e\n',norm(fyyyy - G.y4(x,y), nCh));

H = phi.G(r,s,rx, ry);
fG = H*a;
%fprintf('gradient error = %4.2e\n',norm(fG - G.G(x,y), nCh));

H = phi.L(r,s);
fL = H*a;
%fprintf('Laplacian error = %4.2e\n',norm(fL- G.L(x,y), nCh));

H = phi.B(r,s,rx, ry);
fB = H*a;
%fprintf('biharmonic error = %4.2e\n',norm(fB - G.B(x,y), nCh));

H = phi.D12(r,s,rx,ry);
f12 = H*a;
%fprintf('dx1y2 error = %4.2e\n',norm(f12- G.p12(x,y), nCh));

H = phi.D12(r,s,ry,rx);
f21 = H*a;
%fprintf('dx2y1 error = %4.2e\n',norm(f21- G.p21(x,y), nCh));

H = phi.D22(r,s,ry,rx);
f22 = H*a;
%fprintf('dx2y2 error = %4.2e\n',norm(f22- G.p22(x,y), nCh));

% now interpolate the derivatives at the required point
m=2;
nd=length(cloud);
xy_max=max(max(cloud));
xy_min=min(min(cloud));
r0=5*(xy_max-xy_min)/nd; %basis function radius
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, unkno' );
%RBF_INTERP_ND - does the RBF interpolation
fi = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 1, point' ); %interpolated solution
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, fx );
%RBF_INTERP_ND - does the RBF interpolation
dx = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 1, point' ); %first x-derivative
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, fxx );
%RBF_INTERP_ND - does the RBF interpolation
ddx = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 1, point' ); %second x-derivative
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, fy );
%RBF_INTERP_ND - does the RBF interpolation
dy = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 1, point' ); %first y-derivative
%RBF_WEIGHT - computes the rbf weight values
w = rbf_weight ( m, nd, cloud, r0, @phi1, fyy );
%RBF_INTERP_ND - does the RBF interpolation
ddy = rbf_interp_nd ( m, nd, cloud, r0, @phi1, w, 1, point' ); %second x-derivative