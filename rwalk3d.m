function xyz = rwalk3d(sigma,sz,dim)
% rwalk3d simulates a continuous random walk with a displacement between
% steps of "sigma".

rng('shuffle')
% % Compute the individual steps.
% sz=sz+1;
% xyz=zeros(dim,sz);
% 
% % Stepsize is normal.
% s=sigma*randn(1,sz-1);
% 
% % Direction is random.
% a=randn(dim,sz-1);
% v=s./sqrt(sum(a.^2));
% b=spdiags(v',0,sz-1,sz-1);
% dx(1:dim,1:sz-1)=a*b;
% 
% % Each position is the sum of the previous steps.
% xyz(1:dim,2:sz)=cumsum(dx(1:dim,1:sz-1),2);
% xyz=xyz';

xyz = zeros(sz,dim); % Compute the individual steps
s = sigma*randn(sz,1); % Stepsize is normal
a = randn(sz,dim); % Direction is random
v = s./sqrt(sum(a.^2,2));
b = spdiags(v,0,sz,sz);
dx(1:sz,1:dim) = b*a;
xyz(1:sz,1:dim) = cumsum(dx(1:sz,1:dim),1); % Each position is the sum of the previous steps