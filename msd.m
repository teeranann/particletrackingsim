function xout = msd(x,ran)
% MSD -- function for calculating the mean square displacement in 1D:

% If "x" is an array of positions at different times, then
% dx(t,tau) = x(t+tau) - x(t);
% and
% then <dx> = ave(dx)
% where the ave is over all t values, weighted by the number of values in
% each sum, and 
% <dx^2> = ave(dx^2)
% for plain diffusion,
% <dx^2> = <dy^2> = 2 D tau

% note that the average has to include the fact that the number of
% data points averaged over depends on both t and tau.

% probably not meaningful for tau>length(x)/5 (?)

% no "zero" indices in matlab, so xout(1) is for tau=0
%xlen=length(x);
%ran=1:floor(xlen/5);
%ran=1:xlen;
%ran=1:150;
%tauarray=ran-1;
%ran=1:200;
xout=zeros([size(x,2) length(ran)]);

for cnt=ran,
  tau=cnt+1;
  xttau=x(tau:end,1:end);
  xt=x(1:length(xttau),1:end);
  dxsq=(xttau-xt).^2;
  xout(1:size(x,2),cnt-min(ran)+1)=mean(dxsq,1)';
end