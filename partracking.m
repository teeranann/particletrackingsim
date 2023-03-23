% Particle tracking using the generated fluorescence images
% Teeranan Nongnual
% Random-walk simulations of a particle.
% Generate fluorescence images from particle microsteps. 
% Motion blur could be detected at low shutter speed and fast particle.
% Calculate the deviation (sigma) of the tracked position. 

clear
close all
printtxt=1;
plotflimage=0; %plot 1st and last microstep images to be merged to one motion blur image
plotfluorimg=1; %plot img #1 #2 #3 #last of fluor image and summed
motionblurrwalk=1;
%%
mbstatdyn=1; %static=0 dynamic=1
emrate=150000;
difcoin=5;
motionbh=[];
%%
microstep=10;
simframes=200;
nsamp=1;
%%
for sigmaf=0.100
    %%%%%%%%%
    rmsdsimnm=500;%nm
    rmsdsim=rmsdsimnm/1000;
    te=rmsdsim.^2/6/difcoin;
    for dfrate=1./te %detector framerate in Hz 
        %%
        mbsigma={};
        mbsigma.rmsdsim=[];
        mbsigma.sigmaf=[];
        mbsigma.xloc=[];
        mbsigma.yloc=[];
        mbsigma.sigmaim=[];
        mbsigma.sigmamo=[];
        mbsigma.zxy=[];
        mbsigma.zxybl=[];
        numrmsd=1;
        for sampn=1:nsamp
            fprintf('\n%1.0f/%1.0f\n',sampn,nsamp)
            memallooversz=1;
            while memallooversz==1
                tic
                rmsdsimcur=sqrt(6*difcoin/dfrate)*1000;
                motionb
                if memallooversz==1
                    continue
                end
                trksigma
                if memallooversz==1
                    continue
                end
                if plotfluorimg==1 %summed image
                    figure
                    imgsum=log(sum(imgcnt(:,:,:),3));
                    imagesc(xpsf,ypsf,imgsum); axis equal; colormap jet; hold all;
                    plot(XYZwalksim(:,1),XYZwalksim(:,2),'-w.');hold all
                    plot(Xum(:)+xbdmin-pxlsize/1000,Yum(:)+ybdmin-pxlsize/1000,'-r+');hold all
                    xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
                    axis xy
                    title('summed')
                    ylabel(['y (',char(181),'m)'])
                    xlabel(['x (',char(181),'m)'])
                    background=median(reshape(imgsum(:,:),1,[]));
                    minimg=background;
                    maximg=double(max(max(imgsum)));
                    maximg=ceil(maximg/10^(floor(log(maximg)/log(10))-1))*10^(floor(log(maximg)/log(10)-1));%maximg=470;
                    minimg=floor(minimg);maximg=ceil(maximg);
                    clb=colorbar;
                    set(clb,'Limits',[minimg maximg])
                    set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
                    clb.Label.String='log Counts';
                    clb.Label.FontSize=11;
                    drawnow
                    hold off
                    
                    figure
                    imagesc(xpsf,ypsf,log(sum(imgcnt(:,:,:),3))); axis equal; colormap jet; hold all;
                    xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
                    axis xy
                    %                     title('summed')
                    ylabel(['y (',char(181),'m)'])
                    xlabel(['x (',char(181),'m)'])
                    clb=colorbar;
                    set(clb,'YTick',(minimg:(maximg-minimg)/5:maximg))
                    set(clb,'Limits',[minimg maximg])
                    clb.Label.String='log Counts';
                    set(gca,'fontsize',14,'linewidth',1.2)
                    clb.Label.FontSize=10;
                    clb.FontSize=10;
  
                    set(gca,'tickdir','out')
                    set(gca, 'XTick',min(xpsf):5:max(xpsf),'YTick',min(ypsf):5:max(ypsf))
                    set(gca,'XMinorTick','on','YMinorTick','on')
                    xAx = get(gca,'XAxis');
                    xAx.MinorTickValues=min(xpsf):1:max(xpsf);
                    yAx = get(gca,'yAxis');
                    yAx.MinorTickValues=min(ypsf):1:max(ypsf);
                    ylabel('y')
                    xlabel('x')
                    set(gca,'xticklabel',[],'yticklabel',[])
                    drawnow
                    hold off
                    
                    figure
                    ax11=imagesc(xpsf,ypsf,zeros(size(imgcnt,1), size(imgcnt,2))); axis equal; colormap gray; hold on
                    set(ax11,'visible','off')
                    plot(XYZwalksim(:,1),XYZwalksim(:,2),'-b.');hold all
                    plot(Xum(:)+xbdmin-pxlsize/1000,Yum(:)+ybdmin-pxlsize/1000,'-r.');hold all
                    xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
                    axis xy
                    set(gca,'fontsize',14,'linewidth',1.2)
                    legend('Microsteps','Tracking','fontsize',10,'linewidth',0.8)
                    set(gca,'tickdir','out')
                    set(gca, 'XTick',min(xpsf):5:max(xpsf),'YTick',min(ypsf):5:max(ypsf))
                    set(gca,'XMinorTick','on','YMinorTick','on')
                    xAx = get(gca,'XAxis');
                    xAx.MinorTickValues=min(xpsf):1:max(xpsf);
                    yAx = get(gca,'yAxis');
                    yAx.MinorTickValues=min(ypsf):1:max(ypsf);
                    ylabel('y')
                    xlabel('x')
                    set(gca,'xticklabel',[],'yticklabel',[])
                    drawnow
                    hold off
                    
                    figure;
                    plot(Xwalk,Ywalk,'k','linewidth',0.5);hold on
                    axis equal
                    ylabel('y')
                    xlabel('x')
                    xlim([min(xpsf) max(xpsf)]);ylim([min(ypsf) max(ypsf)])
                    set(gca,'tickdir','out')
                    set(gca, 'XTick',min(xpsf):5:max(xpsf),'YTick',min(ypsf):5:max(ypsf))
                    set(gca,'XMinorTick','on','YMinorTick','on')
                    xAx = get(gca,'XAxis');
                    xAx.MinorTickValues=min(xpsf):1:max(xpsf);
                    yAx = get(gca,'yAxis');
                    yAx.MinorTickValues=min(ypsf):1:max(ypsf);
                    set(gca,'xticklabel',[],'yticklabel',[])
                    set(gca,'fontsize',14,'linewidth',1.2)
                    scatter(Xwalk(1),Ywalk(1),30,'r','filled')
                    scatter(Xwalk(end),Ywalk(end),30,'b','filled')
                    legend('Random-walk','Start','End','fontsize',10,'linewidth',0.8)
                    
                end
                mbsigma.xloc(:,numrmsd)=Xum+xbdmin-pxlsize/1000;
                mbsigma.yloc(:,numrmsd)=Yum+ybdmin-pxlsize/1000;
                if mod(sampn,10)==0
                    fprintf('%1.0i ',sampn)
                    if mod(sampn,100)==0
                        fprintf('\n')
                    end
                end
                sigmaim=sigmaf/sqrt(emrate);
                
                sigmamo=100/sqrt(100000)*(2*rmsdsimcur/sigmaf/1000 - 1).^2;
                fprintf('Sigma theo:\n sigmaim=%1.2fnm, sigmamo=%1.2fnm\n',sigmaim*1000,sigmamo)
                fprintf('Z theo:\n Zxy=%1.6fum^2\n',-4/18*(rmsdsimcur/1000)^2+4*(sigmamo/1000)^2)
                
                mbsigma.rmsdsim(numrmsd)=crmsd*1000; %0.12/sqrt(5000)*sqrt(1+crmsd^2/4/0.12^2)
                mbsigma.sigmaim(numrmsd)=sigmaim*1000;
                mbsigma.sigmamo(numrmsd)=sigmamo;
                mbsigma.sigmaf(numrmsd)=sigmaf*1000;
                
                %%
                %MSD
                XYZwalktblur=[Xum Yum];
                ranmsdosecsim=(1:FA)/dfrate;
                XYZmsdt=msd(XYZwalktblur,1:FA);
                Xmsd=XYZmsdt(1,:);
                Ymsd=XYZmsdt(2,:);
                XYmsdblur=XYZmsdt(1,:)+XYZmsdt(2,:);
                
                [fitX,fitXcoef] = fit(ranmsdosecsim',Xmsd','2*D*x+Z','startpoint',[0.05,0]);
                [fitY,fitYcoef] = fit(ranmsdosecsim',Ymsd','2*D*x+Z','startpoint',[0.05,0]);
                [fitXY,fitXYcoef] = fit(ranmsdosecsim',XYmsdblur','4*D*x+Z','startpoint',[0.05,0]);
                warning('off')
                fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
                fitYv = struct(fitY); fitYcoefv = struct(fitYcoef);
                fitXYv = struct(fitXY); fitXYcoefv = struct(fitXYcoef);
                warning('on')
                [Dx,Zx]=fitXv.coeffValues{1,1:2}; R2x = fitXcoefv.rsquare;
                [Dy,Zy]=fitYv.coeffValues{1,1:2}; R2y = fitYcoefv.rsquare;
                [Dxytrk,Zxytrk]=fitXYv.coeffValues{1,1:2}; R2xy = fitXYcoefv.rsquare;
                
                disp('MSD-TRACKING')
                fprintf(' Dxy=%1.3fum^2/s, Zxy=%1.6fum^2, R2=%1.4f\n',Dxytrk,Zxytrk,R2xy)
                fprintf('Zxytrk-Zxyimg=%1.6fum^2\n',Zxytrk-Zxyimg)
                figure(20)
                plot(ranmsdosecsim',XYmsdsim,'b'); hold all
                plot(fitXYsim,':r');
                plot(ranmsdosecsim',XYmsdblur,'r');
                plot(fitXY,':r')
                mbsigma.zxy(numrmsd)=Zxytrk;
                mbsigma.zxybl(numrmsd)=Zxyimg;
                
                numrmsd=numrmsd+1;
            end
        end
        savename=['mbmsd-s' num2str(sigmaf*1000) '-rmsd' num2str(rmsdsimcur) '.mat'];
        save(savename,'mbsigma');
        fprintf(2,['\nThe output was saved as ' savename '\n'])
    end
    
end