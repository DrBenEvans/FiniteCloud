function CirclePointGenerator
clear
clc
clf
%filename
filename='circle'
%set value of geometric parameters
a=0.0;  %x-coordinate of circle centre
b=0.0;  %y-coordiante of circle centre
c=2.0;  %circle radius
%create the circle points (on  2 curves of 20 points each)
c1=linspace(-pi,0,20);
c2=linspace(pi,0,20);
x1=c*sin(c1);
y1=-c*cos(c1);
x2=c*sin(c2);
y2=c*cos(c2);
x1=x1+a;
x2=x2+a;
y1=y1+b;
y2=y2+b;
scatter(x1,y1)
axis equal
hold on
scatter(x2,y2)
title('Circle Shape')
%write .prf file (.prf)
append1='.prf'
prf_filename=[filename,'',append1];
fileID = fopen(prf_filename,'w');
formatSpec = '%7.4f  %7.4f \n';
for i=1:20
    fprintf(fileID,formatSpec,x1(i),y1(i));
end
for i=2:19
    fprintf(fileID,formatSpec,x2(i),y2(i));
end
fclose(fileID)




