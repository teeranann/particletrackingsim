% Tracking particles in solution and calculate diffusion coefficient by
% 1. MSD
% 2. VACF
% Select particles by using left click. Right click to stop.
% close all
% clear
%%

framesta=1; % Frame number start
frameend=0; % Frame number end % Set '0' to maximum
ranmsdf=20; % Number of frame for plotting VACF and MSD and fitting with diffusion equation

gpxl=ceil((crmsd+fwhmf)/meshsz);
% gpxl=5; % number of pixels (on both sides of center) to use in gaussian fitting
pxl=gpxl+5; % number of pixels (on both sides of center) to find maximum
showtrk=1; % Show tracking video. 1=yes, 0=no.
showevfr=1; %Show picture every n frames. Default is 0 (30 pics).
savevid=1; % Save the tracking video to .avi file. 1=yes, 0=no. %showtrk must be 1 to save vid
FA=20; %Frames to be analysed in MSD and VACF
bp1=1;bp2=30;

%%
pxlsize = meshsz*1000; %65.85nm/pixel

%%
data.nframes=size(imgcnt,3);
data.image=imgcnt;
data.rnddwelltime_s=1/dfrate;

[dim_y,dim_x,numframes] = size(data.image);
background=median(reshape(data.image(:,:,1),1,[]));

bgnoisestd=std(reshape(data.image(:,:,1),1,[]));
% fprintf('Cropped size: %2ix%2ix%2i\n',dim_x,dim_y,numframes)

%%
% figure('units','normalized','position',[0 0 0.50 0.75])
% 
% imagesc(data.image(:,:,1))
minimg=background;
maximg=double(max(max(data.image(:,:,1))));
maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));%maximg=470;
% caxis([minimg maximg])

% axis xy
% set(gca,'dataaspectratio',[1 1 1]) % square pixels
% colormap hot
% clb=colorbar;
% set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
% clb.Label.String='Counts';
% clb.Label.FontSize=11;
% set(gca, 'XTick', (0:2:size(data.image,2)*pxlsize/1000)/pxlsize*1000+1, 'XTickLabel', 0:2:size(data.image,2)*pxlsize/1000)
% set(gca, 'YTick', (0:2:size(data.image,1)*pxlsize/1000)/pxlsize*1000+1, 'YTickLabel', 0:2:size(data.image,1)*pxlsize/1000)
% ylabel(['y (',char(181),'m)'])
% xlabel(['x (',char(181),'m)'])
% drawnow
% hold all

Xg=[];Yg=[];Wg=[];Ag=[];R2g=[];RMSEg=[];Ibcnt=[];Xm=[];Ym=[];

bpimg=fpass(data.image(:,:,1),mkffilt(data.image(:,:,1),bp1,bp2));
[~,maxx]=max(max(bpimg,[],1));
[~,maxy]=max(max(bpimg,[],2));
xidxm1=maxx;
yidxm1=maxy;
Xm(1,end+1)=maxx;
Ym(1,end+1)=maxy;

if yidxm1-pxl<1 || xidxm1-pxl<1 || yidxm1+pxl>size(data.image,1) || xidxm1+pxl>size(data.image,2)
    memallooversz=1;
    disp('Over boundary!')
    return
else
    memallooversz=0;
end
backgroundpar=double(1/4*(data.image(yidxm1-pxl,xidxm1-pxl)+data.image(yidxm1-pxl,xidxm1+pxl)+data.image(yidxm1+pxl,xidxm1-pxl)+data.image(yidxm1+pxl,xidxm1+pxl)));
Ibcnt(1,end+1)=double(sum(sum(data.image((yidxm1-pxl):(yidxm1+pxl),(xidxm1-pxl):(xidxm1+pxl),1))))-(2*pxl+1)^2*backgroundpar;

img=double(data.image(:,:,1));
ctt=cntr2dg(img-background,[xidxm1,yidxm1],gpxl,400/2.35,Ibcnt(1,end)/(2*pi)/(200/65.85)^2); % only do cntr2dgx on raw img
Xg(1,end+1)=ctt.xywi(1,1);
Yg(1,end+1)=ctt.xywi(1,2);
Wg(1,end+1)=ctt.xywi(1,3);
Ag(1,end+1)=ctt.xywi(1,4);
R2g(1,end+1)=ctt.xywi(1,5);
RMSEg(1,end+1)=ctt.xywi(1,6);

% Xg(1,end)=xidxm1;
% Yg(1,end)=yidxm1;

% rectangle('position',[Xg(1,end)-pxl,Yg(1,end)-pxl,2*pxl,2*pxl],'edgecolor','w')

numpar=size(Xg,2);
for par=1:numpar
    for frame=2:numframes
        img=double(data.image(:,:,frame));
        %Maximum
        xidxm_1=round(Xg(frame-1,par));yidxm_1=round(Yg(frame-1,par));

        [~,maxx]=max(max(data.image(yidxm_1-pxl:yidxm_1+pxl,xidxm_1-pxl:xidxm_1+pxl,frame),[],1));
        [~,maxy]=max(max(data.image(yidxm_1-pxl:yidxm_1+pxl,xidxm_1-pxl:xidxm_1+pxl,frame),[],2));
        xidxm=xidxm_1-pxl-1+maxx;
        yidxm=yidxm_1-pxl-1+maxy;

        Xm(frame,par)=xidxm;
        Ym(frame,par)=yidxm;
%         Ibcnt(frame,par)=double(sum(sum(data.image((yidxm-pxl):(yidxm+pxl),(xidxm-pxl):(xidxm+pxl),frame))))-(2*pxl+1)^2*double(background);

        %Gaussian
        ctt=cntr2dg(img-background,[xidxm,yidxm],gpxl,Wg(frame-1,par)/2.35,Ag(frame-1,par)); % Raw img
        Xg(frame,par)=ctt.xywi(1,1);
        Yg(frame,par)=ctt.xywi(1,2);
        Wg(frame,par)=ctt.xywi(1,3);
        Ag(frame,par)=ctt.xywi(1,4);
        R2g(frame,par)=ctt.xywi(1,5);
        RMSEg(frame,par)=ctt.xywi(1,6);
        
%         %maximum
%         Xg(frame,par)=xidxm;
%         Yg(frame,par)=yidxm;
    end
end
% end
% Ib=Ibcnt*camgain/sCMOS_quanteff/collect_eff;
% Idetect=Ibcnt*camgain/sCMOS_quanteff;
% Iphelect=Ibcnt*camgain;

% Show tracking
if showtrk==1
    figure
    if savevid==1
        writerObj = VideoWriter('videoout.avi'); % Name it.
        writerObj.FrameRate = 5; % How many frames per second.
        open(writerObj);
    end
    clr=colormap(lines);
    if showevfr==0
        showevfr=round(numframes/10);
    end
    for frame=1:showevfr:numframes
        if ~isempty('gplot') && ~isempty('maxplot')
            h = findobj('type','line');
            delete(h);
            g = findobj('type','text');
            delete(g);
        end
        imagesc(xpsf,ypsf,data.image(:,:,frame))
        axis xy
        for par=1:numpar
            if showtrk==1
%                 maxplot=scatter(Xm(frame,par),Ym(frame,par),300,'MarkerEdgeColor','r');hold on % Maximum
%                 gplot=scatter(Xg(frame,par),Yg(frame,par),300,'MarkerEdgeColor','w'); hold on % Gaussian
                txfn=text(min(xpsf)+1,min(ypsf)+1.2,['Time = ',sprintf('%1.3f',frame*data.rnddwelltime_s),' s']);
                set(txfn,'color','w');
            else
                gplot=plot(Xg(1:frame,par),Yg(1:frame,par),'color',clr(par,:)); hold on % Gaussian
                maxplot=plot(Xm(1:frame,par),Ym(1:frame,par),'color',clr(2*par,:));
                hold on % Maximum
            end
        end
        caxis([minimg maximg])
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
        colormap jet
        clb=colorbar;
        set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
        clb.Label.String='Counts';
        clb.Label.FontSize=11;
        ylabel(['y (',char(181),'m)'])
        xlabel(['x (',char(181),'m)'])
        
%         imagesc(data.image(:,:,frame))
%         Yax=(1:size(data.image,1))*pxlsize/1000;
%         Xax=(1:size(data.image,2))*pxlsize/1000;
%         set(gca, 'XTick', (0:2:size(data.image,2)*pxlsize/1000)/pxlsize*1000+1, 'XTickLabel', 0:2:size(data.image,2)*pxlsize/1000)
%         set(gca, 'YTick', (0:2:size(data.image,1)*pxlsize/1000)/pxlsize*1000+1, 'YTickLabel', 0:2:size(data.image,1)*pxlsize/1000)
%         xlim([0 max((0:2:size(data.image,2)*pxlsize/1000)/pxlsize*1000+1)])
%         ylim([0 max((0:2:size(data.image,1)*pxlsize/1000)/pxlsize*1000+1)])

        xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
        
        drawnow
        if savevid==1
            fr = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, fr);
        end
    end
    if savevid==1
        close(writerObj);
    end
elseif showtrk==0
%     plot(Xg,Yg)
end

Xum=Xg*pxlsize/1000;
Yum=Yg*pxlsize/1000;
% Wum=Wg*pxlsize/1000;
Zum=Zwalkd(1:size(Xum,1),1);

ranmsdos=ranmsdf*data.rnddwelltime_s;
ranmsdo=round(1:ranmsdos/data.rnddwelltime_s);
ranmsdp=round(1:numframes);
ranmsdosec=ranmsdo*data.rnddwelltime_s;
ranmsdpsec=ranmsdp*data.rnddwelltime_s;
ranmsdpsecsz=size(ranmsdpsec,2);

% for par=1:numpar %Make particles start from (0,0)
%     Xum(:,par)=Xum(:,par)-Xum(1,par);
%     Yum(:,par)=Yum(:,par)-Yum(1,par);
%     %Wum(:,par)=Wum(:,par)-Wum(1,par);
% end

% % MSD
% disp('MSD')
% Xmsd=[];Ymsd=[];Zmsd=[];
% for par=1:numpar
%     Xmsd(:,par)=msd(Xum(ranmsdp,par),ranmsdo);
%     Ymsd(:,par)=msd(Yum(ranmsdp,par),ranmsdo);
%     %     Zmsd(:,par)=msd(Zum(ranmsdp,par),ranmsdo);
% end
% % figure(3)
% Rmsd=[];
% for par=1:numpar
%     Rmsd(:,par)=(Xmsd(:,par)+Ymsd(:,par));
%     [fitX,fitXcoef] = fit(ranmsdosec',Xmsd(:,par),'2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
%     [fitY,fitYcoef] = fit(ranmsdosec',Ymsd(:,par),'2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
%     [fitR,fitRcoef] = fit(ranmsdosec',Rmsd(:,par),'4*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
%     %     [fitZ,fitZcoef] = fit(ranmsdosec',Zmsd(:,par),'2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
%     warning('off')
%     fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
%     fitYv = struct(fitY); fitYcoefv = struct(fitYcoef);
%     fitRv = struct(fitR); fitRcoefv = struct(fitRcoef);
%     %     fitZv = struct(fitZ); fitZcoefv = struct(fitZcoef);
%     warning('on')
%     Dx = fitXv.coeffValues{1,1}; Vx = fitXv.coeffValues{1,2}; Zx = fitXv.coeffValues{1,3}; R2x = fitXcoefv.rsquare;
%     Dy = fitYv.coeffValues{1,1}; Vy = fitYv.coeffValues{1,2}; Zy = fitYv.coeffValues{1,3}; R2y = fitYcoefv.rsquare;
%     Dr = fitRv.coeffValues{1,1}; Vr = fitRv.coeffValues{1,2}; Zr = fitRv.coeffValues{1,3}; R2r = fitRcoefv.rsquare;
% 
% %     plot(ranmsdosec,[Xmsd(:,par),Ymsd(:,par),Rmsd(:,par)]); hold on
% %     plot(fitX,':r')
% %     plot(fitY,':r')
% %     plot(fitR,':r')
% %     grid on
% %     ylabel(['MSDxy (',char(181),'m^2)'])
% %     xlabel('Lag time (s)')
% %     xlim([0 max(ranmsdosec)])
% %     axis square
% 
%     disp([' Particle #' num2str(par)])
%     intrcptx=Zx*1000;
%     intrcpty=Zy*1000;
%     intrcptxy=Zr*1000;
% 
%     fprintf(' Dx=%1.3fum^2/s, V=%1.2eum/s, MSD0x=%1.2fnm, R2=%1.2f\n',Dx,abs(Vx),intrcptx,R2x)
%     fprintf(' Dy=%1.3fum^2/s, V=%1.2eum/s, MSD0y=%1.2fnm, R2=%1.2f\n',Dy,abs(Vy),intrcpty,R2y)
%     fprintf(' Dxy=%1.3fum^2/s, V=%1.2eum/s, MSD0xy=%1.2fnm, R2=%1.2f\n',Dr,abs(Vr),intrcptxy,R2r)
% end
% hold off
% difco=[];
% for par=1:numpar
%     difco(end+1,:)=real([Dx,Dy,Dr,abs(Vx),abs(Vy),abs(Vr),intrcptx,intrcpty,intrcptxy,Zx,Zy,Zr,R2x,R2y,R2r]);
% end
% 
% %VACF
% disp('VACF')
% VACFxt=[];
% VACFyt=[];
% VACFxyt=[];
% % VACFzt=[];
% DvacfxyvsT2t=[];
% DvacfzvsT2t=[];
% for par=1:numpar
%     disp([' Particle #' num2str(par)])
%     XX=Xum(:,par);
%     YY=Yum(:,par);
% %     WW=Wum(:,par);
%     ZZ=Zum(:,par);
%     
%     XYZwalkt=[XX YY ZZ];
%     frate=1/data.rnddwelltime_s;
%     
%     vacfdifco
%     VACFxyt(:,par)=VACFxy;
%     VACFxt(:,par)=VACFx;
%     VACFyt(:,par)=VACFy;
%     %     VACFzt(:,par)=VACFz;
%     DvacfxyvsT2t(:,par)=DvacfxyvsT2;
%     DvacfzvsT2t(:,par)=DvacfzvsT2;
% end
% 
% % sinmul=questdlg('Save trk.mat file and figure?','Save?','No','Yes','No');
% sinmul='No';
% 
% % Output
% if strcmp(sinmul,'Yes')
%     save([fnam2 '-trk.mat'],'Xum','Yum','Xg','Yg','Ag','Wg','Ib','Ibcnt','Wum','Xm','Ym','difco','pxl','gpxl','ranmsdos','framesta','frameend','ranmsdpsec');
%     save([fnam2 '-vacf.mat'],'T1','T2','VACFxyt','VACFxt','VACFyt','VACFzt','DvacfxyvsT2t','DvacfzvsT2t');
%     figure(1)
%     imagesc(data.image(:,:,1))
%     minimg=100;
%     maximg=double(max(max(data.image(:,:,1))));
%     maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));
%     caxis([minimg maximg])
%     axis xy
%     set(gca,'dataaspectratio',[1 1 1]) % square pixels
%     colormap hot
%     clb=colorbar;
%     set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
%     drawnow
%     hold all
%     numpar=size(Xg,2);
%     for par=1:numpar
%         xidx=round(Xg(1,par));
%         yidx=round(Yg(1,par));
%         figure(1)
%         rectangle('position',[xidx-pxl,yidx-pxl,2*pxl,2*pxl],'edgecolor','w')
%         txth=text(xidx+pxl+2,yidx,num2str(par));
%         set(txth,'color',[0 1 0]); % color green
%     end
%     saveas(figure(1),[fnam2 '-img.fig']);
%     saveas(figure(3),[fnam2 '-msdxy.fig']);
%     saveas(figure(4),[fnam2 '-msdz.fig']);
%     saveas(figure(5),[fnam2 '-vacf.fig']);
%     saveas(figure(2),[fnam2 '-xyziw.fig']);
% end
% cd(cdpwd)

% figure(10)
% sumimage=sum(data.image,3);
% imagesc(log(sumimage))
% minimg=log(background*cframes);
% maximg=log(1/2*( max(reshape(sumimage,[],1)) + median(reshape(sumimage,[],1)) ));
% % maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));
% caxis([minimg maximg])
% axis xy
% set(gca,'dataaspectratio',[1 1 1]) % square pixels
% colormap hot
% clb=colorbar;
% set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
% drawnow
% hold all
% numpar=size(Xg,2);
% for par=1:numpar
%     plot(Xg,Yg,'g')
% end