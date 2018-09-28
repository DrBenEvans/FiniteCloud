function [ z ] = Levy( n, m, beta )
%Levy - performs a Levy flight calculation for n steps in m dimensions with
%Power law index Beta
num = gamma(1+beta)*sin(pi*beta/2);
den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma_u = (num/den)^(1/beta);
u = random('Normal',0,sigma_u^2,n,m);
v = random('Normal',0,1,n,m);
z = u./(abs(v).^(1/beta));
end

