function unkno=TIMES(residual,np,nb,unkno,ptype,dt);
%TIMES performs the solution update based on the residual calculation
for ip=nb+1:np
    if(ptype(ip)==3)  %i.e. a domain point
%      dunk=-dt*residual(ip);
      dunk=dt*residual(ip);
      unkno(ip)=unkno(ip)+dunk;
    end
end

