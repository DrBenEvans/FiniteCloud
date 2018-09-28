function [cloud,np,nb,no,ptype,pinsert] = CHDUP(cloud,np,nb,no,ptype,pinsert,tolp)
%CHDUP - checks for duplicate points and removes i.e. points that are
%closer together than the user defined 'tolp' parameter
flag=zeros(1,np);  %zero the point removal flag vector
for ip=1:np
    x=cloud(1,ip);
    y=cloud(2,ip);
    for ich=1:np
        if((ip~=ich) & (flag(ip)==0))
            xch=cloud(1,ich);
            ych=cloud(2,ich);
            dist=sqrt((x-xch)^2+(y-ych)^2);
            if(dist<tolp)  %flag this as point to be removed
                flag(ich)=1;
                    if(ich<=nb)
                      nb=nb-1;
                    elseif(ich<=(nb+no))
                      no=no-1;
                    end
            end
        end
    end
end
% now remove duplicate points from cloud
ipnew=0;
for ip=1:np
    if(flag(ip)==0)
        ipnew=ipnew+1;
        cloud_new(1,ipnew)=cloud(1,ip);
        cloud_new(2,ipnew)=cloud(2,ip);
        ptype_new(ipnew)=ptype(ip);
        pinsert_new(ipnew)=pinsert(ip);
    end
end
np=ipnew;
cloud=cloud_new;
ptype=ptype_new;
pinsert=pinsert_new;
end

