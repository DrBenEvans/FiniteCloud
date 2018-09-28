function [pinsert]=INSERT(np,nb,residual,pinsert,ptins);
%INSERT - this routine updates pinsert based on the residual at each point
%in the cloud
for ip=nb+1,np
    if(abs(residual(ip))>ptins)
        pinsert(ip)=1;
    else
        pinsert(ip)=0;
    end
end

end

