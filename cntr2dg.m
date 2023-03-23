function centers = cntr2dg(img,pks,sz,sigmainput,ainput)

% fitting of 2D gaussian 

pxlsize = 65.85; %65.85nm/pixel

% Set up X and Y variables
% Be careful not to mix up rows/columns (need to check later)
imgsz=size(img);
X=1:imgsz(2); % let X
Y=1:imgsz(1); % let Y

xi=pks(1); % center of gaussian from "peaks"
yi=pks(2);

% find array indices within Region of Interest (ROI)
XROIidx=((xi-sz)<=X) & (X<=(xi+sz));
YROIidx=((yi-sz)<=Y) & (Y<=(yi+sz));

% use only points within ROI for later calculations
XROI=X(XROIidx);
YROI=Y(YROIidx);
%XimgROI=img(circidx);
imgROI=img(YROI,XROI);

imgv=mat2col(imgROI); % flatten to column vector so can use mldivide

x2=XROI*pxlsize;
y2=YROI*pxlsize;
[xx,yy]=meshgrid(x2,y2);

xi=pks(1)*pxlsize; % center of gaussian from "peaks"
yi=pks(2)*pxlsize;

% %% Pre-fit
% % bandpass parameters
% bp1 = 1;
% bp2 = 30;
% H=mkffilt(imgROI,bp1,bp2);
% % use Jason's fourier filter function
% bpimg=fpass(imgROI,H);
% [xData, yData, zData] = prepareSurfaceData( xx, yy, bpimg);
% 
% % Set up fittype and options.
% ft = fittype( 'a1*exp(-(x-xo)^2/(2*sigma^2)-(y-yo)^2/(2*sigma^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf 3 xi-2 yi-2];
% opts.StartPoint = [0 4 xi yi];
% opts.Upper = [Inf 15 xi+2 yi+2];
% 
% % Fit model to data.
% [fitresult, gof] = fit( [xData, yData], zData, ft, opts );
% warning('off')
% fitXv = struct(fitresult); fitXcoefv = struct(gof);
% warning('on')
% %a = fitXv.coeffValues{1,1}; 
% sigmaguess = fitXv.fCoeffValues{1,2};


%% Fit
[xData, yData, zData] = prepareSurfaceData( xx, yy, imgROI);

% Set up fittype and options.
ft = fittype( 'a1*exp((-(x-xo)^2-(y-yo)^2)/(2*sigma*sigma))', 'independent', {'x', 'y'}, 'dependent', 'z');
%ft = fittype( 'a1*exp(-(x-xo)^2/(2*sigma^2)-(y-yo)^2/(2*sigma^2))', 'independent', {'x', 'y'}, 'dependent', 'z');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
minFWHM=100;
maxFWHM=1000; %nm
opts.Lower = [ainput-200 minFWHM/2.35 xi-200 yi-200];
opts.StartPoint = [ainput sigmainput xi yi];
opts.Upper = [ainput+200 maxFWHM/2.35 xi+200 yi+200];
opts.Tolx= 1e-12;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
warning('off')
fitXv = struct(fitresult); fitXcoefv = struct(gof);
warning('on')
a = fitXv.fCoeffValues{1,1}; 
sigma = fitXv.fCoeffValues{1,2}; 
xo = fitXv.fCoeffValues{1,3}; 
yo = fitXv.fCoeffValues{1,4};
R2 = fitXcoefv.rsquare;
RMSE = fitXcoefv.rmse;

centers.xywi=[xo/pxlsize yo/pxlsize 2*sqrt(2*log(2))*sigma/pxlsize a R2 RMSE];