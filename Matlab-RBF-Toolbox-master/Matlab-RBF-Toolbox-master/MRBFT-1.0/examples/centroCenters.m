% centroCenters.m
%
% Distributes Hammersley points in a complexly shaped domain.  The centers
% denser in near the boundary than in the interior.

f = @(t) 3*nthroot( cos(3*t) + sqrt(4 - (sin(3*t)).^2), 3 );    % domain boundary

small = 0.005;
N = 6000;           % N potiential centers in boundary region
boundaryLayerSize = 0.5;

t =  linspace(0,2*pi,200);          
x = f(t).*cos(t);  y = f(t).*sin(t);

% determine the size of the rectangle needed to cover the domain
A = min(x) - small; B = max(x) + small;
C = min(y) - small; D = max(y) + small;

[xc, yc] = rbfCenters.Hammersley2d(N);
xc = (B - A)*xc + A;             % [0,1] --> [A,B]
yc = (D - C)*yc + C;             % [0,1] --> [C,D]

% ---------- boundary region centers -----------------------

th = atan2(yc,xc);  p = sqrt(xc.^2 + yc.^2);
ro = f(th);                   % outter boundary
ri = ro - boundaryLayerSize;  % inner border of boundary region

xn = zeros(N,1);  yn = zeros(N,1);   I = 1;
for i=1:N
    if  and( p(i) < ro(i), p(i)>ri(i) )
        xn(I) = xc(i); yn(I) = yc(i);  I = I + 1;
    end
end
xc = xn(1:find(xn,1,'last'));    % remove trailing zeros
yc = yn(1:find(yn,1,'last'));

% ------ interior centers----------------------------

N2 = 2100;     % N2 potiential centers in interior region
[xci, yci] = rbfCenters.Hammersley2d(N2);

A = A + boundaryLayerSize;  B = B - boundaryLayerSize;
C = C + boundaryLayerSize;  D = D - boundaryLayerSize;
xci = (B - A)*xci + A;             % [0,1] --> [A,B]
yci = (D - C)*yci + C;             % [0,1] --> [C,D]

th = atan2(yci,xci); p = sqrt(xci.^2 + yci.^2);
ri = f(th) - boundaryLayerSize;    % interior region boundary

xn = zeros(N2,1);  yn = zeros(N2,1);  I = 1;
for i=1:N2
    if  p(i)<(ri(i) - small) 
        xn(I) = xci(i); yn(I) = yci(i);  I = I + 1;
    end
end
 
xci = xn(1:find(xn,1,'last'));    % remove trailing zeros
yci = yn(1:find(yn,1,'last'));
x = [xc; xci];   y = [yc; yci];   % merge centers


% --find centers in half of the domain (x-axis symmetry) ----------

I = find(y>(0 + 0e-3 )); x = x(I);  y = y(I);

% ------ extend "centro-symmetrically" to the other half -----------

x = [x; flipud(x)];   y = [y; flipud(-y)];
 
% ---------------- verify centrosymmetry --------------------------

r = abs(rbfx.distanceMatrix2d(x,y));
centro = rbfCentro.isCentro(r)               

% -----------------------------------------------------------------
scatter(x,y,'b.')