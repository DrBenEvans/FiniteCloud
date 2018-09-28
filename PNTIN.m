function [cloud,np,ptype,pinsert] = PNTIN(cloud,np,nb,no,ptype,pinsert,pnorm,nLevy,beta,bstep)
%PNTIN - this routine uses a Levy flight approach to insert points into the point cloud
for ip=1:np
    np_new=np;
    if(pinsert(ip)==1)
        np_new=np_new+1;
        z = Levy(nLevy,2,beta);
        if((ptype(ip)==1)|(ptype(ip)==2)) %if boundary point
          cloud(1,np_new)=cloud(1,ip) + bstep*pnorm(1,ip); %first step away from boundary in the normal direction
          cloud(2,np_new)=cloud(2,ip) + bstep*pnorm(2,ip);
          for istep=1:nLevy  %perform Levy flight steps (from boundary points) and check whether intersecting the boundary at each step
            [flag, inter] = CHINT([cloud(1,np_new),cloud(2,np_new)],[cloud(1,np_new) + 0.1*z(istep,1),cloud(2,np_new) + 0.1*z(istep,2)],cloud,nb,no);
            if(flag==0)  %no intersect
              cloud(1,np_new)=cloud(1,np_new) + 0.1*z(istep,1);
              cloud(2,np_new)=cloud(2,np_new) + 0.1*z(istep,2);
            end
          end          
        else
          cloud(1,np_new)=cloud(1,ip); %initialise new point
          cloud(2,np_new)=cloud(2,ip);
          for istep=1:nLevy %perform Levy flight steps and check whether intersecting the boundary at each step
            [flag, inter] = CHINT([cloud(1,np_new),cloud(2,np_new)],[cloud(1,np_new) + 0.1*z(istep,1),cloud(2,np_new) + 0.1*z(istep,2)],cloud,nb,no);
            if(flag==0) %no intersect
              cloud(1,np_new)=cloud(1,np_new)+0.1*z(istep,1);
              cloud(2,np_new)=cloud(2,np_new)+0.1*z(istep,2);  %NEED A CHECK AT THIS POINT TO SEE IF INSERTED POINT IS INSIDE THE SHAPE!
            end
          end
        end
        ptype(np_new)=3;  %internal point
        pinsert(np_new)=0;
    end
    np=np_new;
end

