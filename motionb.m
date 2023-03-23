%%
if ~exist('motionbh','var')
    plotflimage=0;
    printtxt=1;
    
    dfrate=1000; %detector framerate in Hz
    microstep=100;
    simframes=200;
    difcoin=5;
    emrate=100000;
    sigmaf=0.100;
    rmsdsimcur=sqrt(6*difcoin/dfrate);
end

%% r-walk parameters

rwalk;
Xwalkd=Xwalk;Ywalkd=Ywalk;Zwalkd=Zwalk; %xyz in um

%% Detector parameters

pxlsz
meshsz=pxlsize/1000; %detector pixel size in um
cfrate=frate/dfrate; %collecting framerate
% cfrate=floor(frate/dfrate);
% dfrate=frate/cfrate;
zdof=30; %depth of field in um
readout=2;%1.7;
gain=0.6;
backgr=100;

% emrate=1000; %emissionrate cnt/dframe
cntppar=emrate/cfrate; %cntppar=1e5/dfrate; %emissionrate/dfrate in cnt/par/frame
% sigmaf=0.100; %gaussian sigma in um of focused psf
cntmaxf=cntppar/(2*pi*sigmaf^2); %Imax=Itotal/(2*pi*sigma^2) (gaussian dist.)
cntmaxppxf=cntmaxf*(meshsz)^2;
fwhmf=sigmaf*2.355;
cframes=floor(nframes/cfrate); %number of collected frames

% if cframes>2000
%     cframes=2000;
%     Xwalkd=Xwalkd-mean(Xwalkd(1:ceil(cframes*cfrate)));
%     Ywalkd=Ywalkd-mean(Ywalkd(1:ceil(cframes*cfrate)));
%     Zwalkd=Zwalkd-mean(Zwalkd(1:ceil(cframes*cfrate))); %xyz in um
% end

crmsd=sqrt(2*dim*difcoin/dfrate);


%%
xsizeum=1*(max(Xwalkd(1:ceil(cframes*cfrate)))-min(Xwalkd(1:ceil(cframes*cfrate))) ) +1;
ysizeum=1*(max(Ywalkd(1:ceil(cframes*cfrate)))-min(Ywalkd(1:ceil(cframes*cfrate))) ) +1;
xbdmax=ceil(max(Xwalkd))+2*ceil(crmsd+fwhmf); xbdmin=floor(min(Xwalkd))-2*ceil(crmsd+fwhmf); %XY boundary 5um away from XYwalk
ybdmax=ceil(max(Ywalkd))+2*ceil(crmsd+fwhmf); ybdmin=floor(min(Ywalkd))-2*ceil(crmsd+fwhmf);
% xbdmax=ceil(max(Xwalkd))+xsizeum+2*rmsdsimcur; xbdmin=floor(min(Xwalkd))-xsizeum-2*rmsdsimcur; %XY boundary 5um away from XYwalk
% ybdmax=ceil(max(Ywalkd))+ysizeum+2*rmsdsimcur; ybdmin=floor(min(Ywalkd))-ysizeum-2*rmsdsimcur;

xpsf=xbdmin:meshsz:xbdmax; ypsf=ybdmin:meshsz:ybdmax; %XY vector 
[xpsfgrid,ypsfgrid]=meshgrid(xpsf,ypsf); %XY grid (x:column,y:row)
imgszallo=size(ypsf,2)*size(xpsf,2)*cframes/1024^3*8;
mempreallo=memory;
if printtxt~=0
    fprintf(2,'!Imagesizeprecal=%1.2fGB/%1.2fGB\n',imgszallo,mempreallo.MemAvailableAllArrays/1024^3)
end
if imgszallo > mempreallo.MemAvailableAllArrays/1024^3
    memallooversz=1;
    disp('Oversized!')
    return
else
    memallooversz=0;
end

imgcnt=(zeros(size(ypsf,2),size(xpsf,2),cframes));
% imgcntpf=zeros(size(ypsf,2),size(xpsf,2),cframes);
dtt=whos('imgcnt');memav=memory;

if printtxt~=0
    fprintf('Imagesize=%1.2fGB/%1.2fGB\n',dtt.bytes*9.53674e-7/1024,memav.MemAvailableAllArrays/1024^3)
    fprintf('Detector parameters:\n F=%1.4fHz, N=%1.2fframes, Time=%1.2fs, RMSD=%1.3fum, L=%1.3fum,',dfrate,cframes,cframes/dfrate,crmsd,sqrt(delta/dfrate*cframes))
end

% %Ne
% for i=1:nframes
%     fwhmz=sqrt(((Zwalkd(i))^2/zdof+1)^2)*(fwhmf)^2;
%     sigmaz=fwhmz*2.355;
%     cntmaxz=cntperpar/(2*pi*sigmaz^2);
%     cntmaxperpxz=cntmaxz*(meshsz/1000)^2;
%     
% %     imgcnt(:,:,i)=zeros(size(xpsf,2),size(ypsf,2));
% %     [~,xmax]=min(abs(Xwalkd(i)-xpsf));[~,ymax]=min(abs(Ywalkd(i)-ypsf));
% %     xspotmax=xmax+2/meshsz*1000;xspotmin=xmax-2/meshsz*1000;
% %     yspotmax=ymax+2/meshsz*1000;yspotmin=ymax-2/meshsz*1000;
% %     imgcnt(xspotmin:xspotmax,yspotmin:yspotmax,i)=poissrnd((cntmaxperpxz*exp((-(xpsfgrid(xspotmin:xspotmax,yspotmin:yspotmax)-Xwalkd(i)).^2-(ypsfgrid(xspotmin:xspotmax,yspotmin:yspotmax)-Ywalkd(i)).^2)/(2*sigmaz^2))));
% 
%     imgcnt(:,:,i)=poissrnd((cntmaxperpxz*exp((-(xpsfgrid-Xwalkd(i)).^2-(ypsfgrid-Ywalkd(i)).^2)/(2*sigmaz^2))));
%     
% %     psf(:,:,i)=sqrt();
% %     psf(:,:,i)=hist3([Xwalkd((i-1)*collectfrate+1:i*collectfrate),Ywalkd((i-1)*collectfrate+1:i*collectfrate)],'edges',{xboundarymin:meshsize/1000:xboundarymax,yboundarymin:meshsize/1000:yboundarymax});
% end
% % %Readout noise
% % for i=1:nframes
% %     fwhmz=sqrt(((Zwalkd(i))^2/zdof+1)^2)*(fwhmf)^2;
% %     sigmaz=fwhmz*2.355;
% %     psf(:,:,i)=poissrnd(1000*exp((-(xpsfgrid-Xwalkd(i)).^2-(ypsfgrid-Ywalkd(i)).^2)/(2*sigmaz^2)));
% % end
% % %Ne+Readout then poisson
% % for i=1:nframes
% %     imgcnt(:,:,i)=poissrnd(psf(:,:,i)+reado(:,:,i));
% % end
% 
% figure
% imagesc(xpsf,ypsf,imgcnt(:,:,1)); axis equal; colormap hot;

% xpsfs=xboundarymin:meshsz/1000:xboundarymax;
% ypsfs=yboundarymin:meshsz/1000:yboundarymax;
imgcntm=zeros(size(ypsf,2),size(xpsf,2),microstep);
for i=1:cframes
    
    microscnt=1;
    for j=round( (i-1)*cfrate+1:i*cfrate )
        fwhmz=sqrt(((Zwalkd(j))^2/zdof+1)^2)*(fwhmf)^2;
        sigmaz=fwhmz*2.355;
        cntmaxz=cntppar/(2*pi*sigmaz^2);
        cntmaxppxz=cntmaxz*(meshsz)^2;
        
%         [~,xmax]=min(abs(Xwalkd(i)-xpsf));[~,ymax]=min(abs(Ywalkd(i)-ypsf));
%         swidth=round(3*(0.25+crmsd)/meshsz);
%         xspotmax=xmax+swidth;xspotmin=xmax-swidth;
%         yspotmax=ymax+swidth;yspotmin=ymax-swidth;
%         imgcntbox(yspotmin:yspotmax,xspotmin:xspotmax,i)=((cntmaxperpxz*exp((-(xpsfgrid(yspotmin:yspotmax,xspotmin:xspotmax)-Xwalkd(i)).^2-(ypsfgrid(yspotmin:yspotmax,xspotmin:xspotmax)-Ywalkd(i)).^2)/(2*sigmaz^2))));
%         imgcnt(yspotmin:yspotmax,xspotmin:xspotmax,i)=imgcnt(yspotmin:yspotmax,xspotmin:xspotmax,i)+imgcntbox(yspotmin:yspotmax,xspotmin:xspotmax,i);
        imgcntm(:,:,microscnt)=((cntmaxppxf*exp((-(xpsfgrid-Xwalkd(j)).^2-(ypsfgrid-Ywalkd(j)).^2)/(2*sigmaz*sigmaz))));
        microscnt=microscnt+1;
%         imgcnt(:,:,i)=imgcnt(:,:,i)+((cntmaxperpxz*exp((-(xpsfgrid-Xwalkd(j)).^2-(ypsfgrid-Ywalkd(j)).^2)/(2*sigmaz*sigmaz))));
        if plotflimage==1
            %PLOT image at frame #1 and #end before summing
            if (j == 1 || j == microstep) && i==1
                figure
                imagesc(xpsf,ypsf,imgcntm(:,:,j)); axis equal; colormap jet; hold all;
%                 set(get(colorbar,'label'),'string','x 10^3')
                axis equal
                axis xy
                xlim([-10 10])
                ylim([-10 10])
                ylabel(['y (',char(181),'m)'])
                xlabel(['x (',char(181),'m)'])

            end
        end
    end
    %no poisson
%     imgcnt(:,:,i)=round(sum(imgcntm,3)+round(((normrnd(0,readout/gain,size(ypsf,2),size(xpsf,2)))+(ones(size(ypsf,2),size(xpsf,2))*backgr)))); % x=column,y=row
    %poisson
    imgcnt(:,:,i)=round(poissrnd(sum((imgcntm),3))+(((normrnd(0,readout/gain,size(ypsf,2),size(xpsf,2)))+(ones(size(ypsf,2),size(xpsf,2))*backgr)))); % x=column,y=row
end
pxl=2*ceil(crmsd/pxlsize*1000)+2;
[~,maxx]=max(max(imgcnt(:,:,1),[],1));
[~,maxy]=max(max(imgcnt(:,:,1),[],2));
cntperparf1=double(sum(sum(imgcnt((maxy-pxl):(maxy+pxl),(maxx-pxl):(maxx+pxl),1))))-(2*pxl+1)^2*backgr;
if printtxt~=0
    fprintf(' max(I)=%1.0fcnt/px/frame, I-b=%1.0fcnt/par/frame\n',max(reshape(imgcnt,1,[])),cntperparf1)
end
if plotfluorimg
    %PLOT fluorescence image (summed all spots in microsteps)
    for k=[1 size(imgcnt,3)]
        figure
        imagesc(xpsf,ypsf,(imgcnt(:,:,k))); axis equal; colormap jet; hold all;
        axis equal
        axis xy
        xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
        % plot(0,0,'r+')
%         plot([-1000 1000],[0 0],'w')
%         plot([0 0],[-1000 1000],'w')
%         xlim([-10 10])
%         ylim([-10 10])
        ylabel(['y (',char(181),'m)'])
        xlabel(['x (',char(181),'m)'])
        title(k)
        background=median(reshape(imgcnt(:,:,1),1,[]));
        minimg=background;
        maximg=double(max(max(imgcnt(:,:,1))));
        maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));%maximg=470;
        clb=colorbar;
        set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
        clb.Label.String='Counts';
        clb.Label.FontSize=11;
        drawnow
    end
end

%% msd fluor image
Xwalkimg=Xwalk(microstep:microstep:microstep*simframes);
Ywalkimg=Ywalk(microstep:microstep:microstep*simframes);
XYZwalksim=[Xwalkimg Ywalkimg];
ranmsdosecsim=(1:FA)/dfrate;
XYZmsdt=msd(XYZwalksim,1:FA);
Xmsd=XYZmsdt(1,:);
Ymsd=XYZmsdt(2,:);
XYmsdsim=XYZmsdt(1,:)+XYZmsdt(2,:);

if mbstatdyn==1 %static=0 dynamic=1
    [fitXsim,fitXcoef] = fit(ranmsdosecsim',Xmsd','2*D*x+Z','startpoint',[0.05,0]);
    [fitYsim,fitYcoef] = fit(ranmsdosecsim',Ymsd','2*D*x+Z','startpoint',[0.05,0]);
    [fitXYsim,fitXYcoef] = fit(ranmsdosecsim',XYmsdsim','4*D*x+Z','startpoint',[0.05,0]);
    warning('off')
    fitXv = struct(fitXsim); fitXcoefv = struct(fitXcoef);
    fitYv = struct(fitYsim); fitYcoefv = struct(fitYcoef);
    fitXYv = struct(fitXYsim); fitXYcoefv = struct(fitXYcoef);
    warning('on')
    [Dx,Zx]=fitXv.coeffValues{1,1:2}; R2x = fitXcoefv.rsquare;
    [Dy,Zy]=fitYv.coeffValues{1,1:2}; R2y = fitYcoefv.rsquare;
    [Dxyimg,Zxyimg]=fitXYv.coeffValues{1,1:2}; R2xy = fitXYcoefv.rsquare;
    
    
%     fprintf(' Dx=%1.3fum^2/s, V=%1.2eum/s, MSD0x=%1.2fnm, R2=%1.2f\n',Dx,abs(Vx),sqrt(abs(Zx)/2)*1000,R2x)
%     fprintf(' Dy=%1.3fum^2/s, V=%1.2eum/s, MSD0y=%1.2fnm, R2=%1.2f\n',Dy,abs(Vy),sqrt(abs(Zy)/2)*1000,R2y)
%     fprintf(' Dz=%1.3fum^2/s, V=%1.2eum/s, MSD0z=%1.2fnm, R2=%1.2f\n',Dz,abs(Vz),sqrt(abs(Zz)/2)*1000,R2z)
%     fprintf(' Dxy=%1.3fum^2/s, V=%1.2eum/s, MSD0xy=%1.2fnm, R2=%1.2f\n',Dxy,abs(Vxy),sqrt(abs(Zxy)/2)*1000,R2xy)
%     fprintf(' Dxyz=%1.3fum^2/s, V=%1.2eum/s, MSD0r=%1.2fnm, R2=%1.2f\n',Dr,abs(Vr),sqrt(abs(Zr)/2)*1000,R2r)
    
    disp('MSD-IMAGE')
    fprintf(' Dxy=%1.3fum^2/s, Zxy=%1.6fum^2, R2=%1.4f\n',Dxyimg,Zxyimg,R2xy)
    figure(20)
    plot(ranmsdosecsim',XYmsdsim); hold all
    %         plot(fitX,':r');
    %         plot(fitY,':r')
    plot(fitXYsim,':r'); hold all
end



%%

% minimg=background;
% maximg=double(max(max(imgcnt(:,:,1))));
% maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));%maximg=470;
% caxis([minimg maximg])
% axis xy
% set(gca,'dataaspectratio',[1 1 1]) % square pixels
% colormap hot
% clb=colorbar;
% set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
% clb.Label.String='Counts';
% clb.Label.FontSize=11;
% set(gca, 'XTick', (0:2:size(imgcnt,2)*pxlsize/1000)/pxlsize*1000+1, 'XTickLabel', 0:2:size(imgcnt,2)*pxlsize/1000)
% set(gca, 'YTick', (0:2:size(imgcnt,1)*pxlsize/1000)/pxlsize*1000+1, 'YTickLabel', 0:2:size(imgcnt,1)*pxlsize/1000)

% 
% h=implay(imgcnt,10);
% h.Visual.setPropertyValue('UseDataRange',true);
% h.Visual.setPropertyValue('DataRangeMin',100);
% h.Visual.setPropertyValue('DataRangeMax',max(max(max(imgcnt))));
% h.Visual.ColorMap.Map=jet(256);
