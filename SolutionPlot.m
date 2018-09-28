%quick and dirty script to plot the solution of the rectangular,
%sin(pi*x/a) problem
for i=1:np
    solution(i)=(1/sinh(pi))*sin(pi*cloud(1,i))*sinh(pi*cloud(2,i));
end
% visualise the residual distribution
   tri=delaunay(cloud(1,:),cloud(2,:));
   figure(100)
   trisurf(tri,cloud(1,:),cloud(2,:),solution,solution)
   title('Delaunay/trisurf representation of the solution')
   view(0,90)
   hold on
   scatter(cloud(1,1:nb),cloud(2,1:nb),'r')
   axis equal
   colorbar