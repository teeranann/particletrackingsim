% walkfit is for fitting Xum Yum from walk2d.m
close all

plotmsd=0; plotMSDVACF=0; plot2D3D=0;
saveXYZwalkTs=0;

% difcoin=2; % um^2/s
frate=dfrate*microstep; % Hz
nframes=frate*simframes/dfrate;
sz=[1 nframes]; % Size of random walk

dim = 3; % dimension (walk2d uses 2D)
delta = 2*dim*difcoin; % MSD/sec
sigma = sqrt(delta/frate); % RMSD/step
FA = 20; % Frame number to be analyzed
ranmsdosec = (1:FA)/frate;
Ts=(1:sz(2))/frate;

confined=0; % Free diff =0, Confined diff =1;
rsphere=0.125; %um
confnode=10; %must be even numbers
confraction=0.99;

% fprintf('RMSD=%1.3fum 6*RMSD=%1.3fum Pore=%1.3fum\n',sigma,6*sigma,2*rsphere)
if printtxt~=0
    fprintf(2,'\nR-Walk INPUT:\n D=%1.3fum^2/s, F=%1.2fHz, N=%1.0fframes, Time=%1.2fs, RMSD=%1.3fum, L=%1.3fum\n',difcoin,frate,sz(2),sz(2)/frate,sigma,sqrt(delta/frate*sz(2)))
end
reflectz=0;

sp=0; % Small pertubating force
difcosp=0.01; % um^2/s
deltasp = 2*dim*difcosp; % MSD/sec
sigmasp = sqrt(deltasp/frate); % RMSD/step

%% RWALK
XYZwalk=rwalk3d(sigma,round(sz(2)),dim); % micron

if sp==1
    XYZwalksp=rwalk3d(sigmasp,sz(2),dim); % micron
    XYZwalk=XYZwalk+XYZwalksp;
end
XYZwalk=[XYZwalk(1:end,1)-XYZwalk(1,1) XYZwalk(1:end,2)-XYZwalk(1,2) XYZwalk(1:end,3)-XYZwalk(1,3)];
Xwalk=XYZwalk(1:end,1); %micron
Ywalk=XYZwalk(1:end,2);
Zwalk=XYZwalk(1:end,3);

if reflectz==1
    ZwalkNorm=Zwalk; % Keep the old Zwalk
    rng('shuffle')
    mirror=min(Zwalk)+(max(Zwalk)-min(Zwalk)).*rand(1,1);
    Zwalk=-abs(Zwalk-mirror)+mirror;
end

if confined==1
    % Node
    intv=round(length(XYZwalk)/confnode);
    FF=1:length(XYZwalk);
    % Boundary
    frameconfined=[];framebounce=[];
    for j=1:1:confnode/2
        confF=FF((j-1)*2*intv+1:((j-1)*2*intv)+2*intv*confraction);
        Xmove=Xwalk(confF(1));Ymove=Ywalk(confF(1));Zmove=Zwalk(confF(1));
        Xwalk=Xwalk-Xmove;
        Ywalk=Ywalk-Ymove;
        Zwalk=Zwalk-Zmove;
        frametotal=2*intv*confraction;
        frameconfined(j)=size(confF,2);
        framebounce(j)=0;
        for i=confF
            
            if Xwalk(i)^2+Ywalk(i)^2+Zwalk(i)^2 > rsphere^2
                framebounce(j)=framebounce(j)+1;
            end
            
            while Xwalk(i)^2+Ywalk(i)^2+Zwalk(i)^2 > rsphere^2
                X0=Xwalk(i-1);Y0=Ywalk(i-1);Z0=Zwalk(i-1);
                X1=Xwalk(i);Y1=Ywalk(i);Z1=Zwalk(i);
                
                SS = sigma*randn(sz(2),1); % Stepsize is normal
                DI = randn(sz(2),dim); % Direction is random
                v = SS./sqrt(sum(DI.^2,2));
                b = spdiags(v,0,sz(2),sz(2));
                dx(1:sz(2),1:dim) = b*DI;
                XYZf(1:sz(2),1:dim) = cumsum(dx(1:sz(2),1:dim),1);
                
                Xwalk(i-1:sz(2),:)=XYZf(i-1:sz(2),1)-XYZf(i-1,1)+X0;
                Ywalk(i-1:sz(2),:)=XYZf(i-1:sz(2),2)-XYZf(i-1,2)+Y0;
                Zwalk(i-1:sz(2),:)=XYZf(i-1:sz(2),3)-XYZf(i-1,3)+Z0;
                
%                 Xwalk(i,:)=XYZf(i,1)-XYZf(i-1,1)+X0;
%                 Ywalk(i,:)=XYZf(i,2)-XYZf(i-1,2)+Y0;
%                 Zwalk(i,:)=XYZf(i,3)-XYZf(i-1,3)+Z0;
%                 Xwalk(i+1:sz(2),:)=Xwalk(i+1:sz(2),:)-Xwalk(i+1,:)+X1;
%                 Ywalk(i+1:sz(2),:)=Ywalk(i+1:sz(2),:)-Ywalk(i+1,:)+Y1;
%                 Zwalk(i+1:sz(2),:)=Zwalk(i+1:sz(2),:)-Zwalk(i+1,:)+Z1;
            end
        end
        Xwalk=Xwalk+Xmove;
        Ywalk=Ywalk+Ymove;
        Zwalk=Zwalk+Zmove;
    end
end

Xwalk=Xwalk-mean(Xwalk);
Ywalk=Ywalk-mean(Ywalk);
Zwalk=Zwalk-mean(Zwalk);
% Xum=Xwalk;
% Yum=Ywalk;
% Zum=Zwalk;
% % trksol_trajectplot

if confined==1
    fprintf('Total frame = %1.0f, #Pore = %1.0f\n',sz(2),confnode/2)
    for m=1:size(framebounce,2)
        fprintf('Pore %1i: #boucingf/#totalf = %1.0f/%1.0f\n',m,framebounce(m),frametotal)
    end
end


%% MSD

lengthav=mean([max(Xwalk)-min(Xwalk),max(Ywalk)-min(Ywalk),max(Zwalk)-min(Zwalk)]);
if printtxt~=0
    fprintf('R-Walk OUTPUT:\n')
    fprintf(' L=%1.3fum\n',lengthav)
end
XYZwalkt=[Xwalk Ywalk Zwalk];
tic
XYZmsdt=msd(XYZwalkt,1:FA);
if exist('msdtime','var')
    msdtime(end+1)=toc;
end
Xmsd=XYZmsdt(1,:);
Ymsd=XYZmsdt(2,:);
Zmsd=XYZmsdt(3,:);
XYmsd=XYZmsdt(1,:)+XYZmsdt(2,:);
XYZmsd=XYZmsdt(1,:)+XYZmsdt(2,:)+XYZmsdt(3,:);

%2D
if plot2D3D==1
    
    figure
    plot((1:sz(2))/frate,[Xwalk,Ywalk,Zwalk])
    
    figure
    plot(Xwalk,Ywalk); hold all
    axis equal
    axis xy
    scatter(Xwalk(1),Ywalk(1))
    axis equal
    axis xy
    scatter(Xwalk(end),Ywalk(end))
    axis equal
    axis xy
    
    
end

if confined==0
    if printtxt~=0
        [fitX,fitXcoef] = fit(ranmsdosec',Xmsd','2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
        [fitY,fitYcoef] = fit(ranmsdosec',Ymsd','2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
        [fitZ,fitZcoef] = fit(ranmsdosec',Zmsd','2*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
        [fitXY,fitXYcoef] = fit(ranmsdosec',XYmsd','4*D*x+Z','startpoint',[0.05,0]);
        [fitR,fitRcoef] = fit(ranmsdosec',XYZmsd','6*D*x+(V*x)^2+Z','startpoint',[0.05,1e-10,0]);
        warning('off')
        fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
        fitYv = struct(fitY); fitYcoefv = struct(fitYcoef);
        fitZv = struct(fitZ); fitZcoefv = struct(fitZcoef);
        fitXYv = struct(fitXY); fitXYcoefv = struct(fitXYcoef);
        fitRv = struct(fitR); fitRcoefv = struct(fitRcoef);
        warning('on')
        [Dx,Vx,Zx]=fitXv.coeffValues{1,1:3}; R2x = fitXcoefv.rsquare;
        [Dy,Vy,Zy]=fitYv.coeffValues{1,1:3}; R2y = fitYcoefv.rsquare;
        [Dz,Vz,Zz]=fitZv.coeffValues{1,1:3}; R2z = fitZcoefv.rsquare;
        [Dxymicro,Zxymicro]=fitXYv.coeffValues{1,1:2}; R2xy = fitXYcoefv.rsquare;
        [Dr,Vr,Zr]=fitRv.coeffValues{1,1:3}; R2r = fitRcoefv.rsquare;
        

%     fprintf(' Dx=%1.3fum^2/s, V=%1.2eum/s, MSD0x=%1.2fnm, R2=%1.2f\n',Dx,abs(Vx),sqrt(abs(Zx)/2)*1000,R2x)
%     fprintf(' Dy=%1.3fum^2/s, V=%1.2eum/s, MSD0y=%1.2fnm, R2=%1.2f\n',Dy,abs(Vy),sqrt(abs(Zy)/2)*1000,R2y)
%     fprintf(' Dz=%1.3fum^2/s, V=%1.2eum/s, MSD0z=%1.2fnm, R2=%1.2f\n',Dz,abs(Vz),sqrt(abs(Zz)/2)*1000,R2z)
%     fprintf(' Dxy=%1.3fum^2/s, V=%1.2eum/s, MSD0xy=%1.2fnm, R2=%1.2f\n',Dxy,abs(Vxy),sqrt(abs(Zxy)/2)*1000,R2xy)
%     fprintf(' Dxyz=%1.3fum^2/s, V=%1.2eum/s, MSD0r=%1.2fnm, R2=%1.2f\n',Dr,abs(Vr),sqrt(abs(Zr)/2)*1000,R2r)
    
        disp('MSD-MICROSTEPS')
%         fprintf(' Dxyz=%1.3fum^2/s, R2=%1.2f\n',Dr,R2r)
        fprintf(' Dxy=%1.3fum^2/s, Zxy=%1.6fum^2, R2=%1.4f\n',Dxymicro,Zxymicro,R2xy)
    end
elseif confined==1
    [fitX,fitXcoef] = fit(ranmsdosec',Xmsd','L*(1-Z*exp(-2/3*D*x/L))','startpoint',[difcoin,rsphere^2,1]);
    [fitY,fitYcoef] = fit(ranmsdosec',Ymsd','L*(1-Z*exp(-2/3*D*x/L))','startpoint',[difcoin,rsphere^2,1]);
    [fitZ,fitZcoef] = fit(ranmsdosec',Zmsd','L*(1-Z*exp(-2/3*D*x/L))','startpoint',[difcoin,rsphere^2,1]);
    [fitXY,fitXYcoef] = fit(ranmsdosec',XYmsd','L*(1-Z*exp(-2*D*x/L))','startpoint',[difcoin,rsphere^2,1]);
    [fitR,fitRcoef] = fit(ranmsdosec',XYZmsd','L*(1-Z*exp(-6/3*D*x/L))','startpoint',[difcoin,rsphere^2,1]);
    warning('off')
    fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
    fitYv = struct(fitY); fitYcoefv = struct(fitYcoef);
    fitZv = struct(fitZ); fitZcoefv = struct(fitZcoef);
    fitXYv = struct(fitXY); fitXYcoefv = struct(fitXYcoef);
    fitRv = struct(fitR); fitRcoefv = struct(fitRcoef);
    warning('on')
    [Dx,Lx,Zx]=fitXv.coeffValues{1,1:3}; R2x = fitXcoefv.rsquare;
    [Dy,Ly,Zy]=fitYv.coeffValues{1,1:3}; R2y = fitYcoefv.rsquare;
    [Dz,Lz,Zz]=fitZv.coeffValues{1,1:3}; R2z = fitZcoefv.rsquare;
    [Dxymicro,Lxy,Zxymicro]=fitXYv.coeffValues{1,1:3}; R2xy = fitXYcoefv.rsquare;
    [Dr,Lr,Zr]=fitRv.coeffValues{1,1:3}; R2r = fitRcoefv.rsquare;
    disp('MSD: Confined D')
    fprintf(' Dx=%1.3fum^2/s, Lx=%1.3fum, Ax=%1.2f, R2=%1.4f\n',Dx,sqrt(Lx*3),Zx,R2x)
    fprintf(' Dy=%1.3fum^2/s, Ly=%1.3fum, Ay=%1.2f, R2=%1.4f\n',Dy,sqrt(Ly*3),Zy,R2y)
    fprintf(' Dz=%1.3fum^2/s, Lz=%1.3fum, Az=%1.2f, R2=%1.4f\n',Dz,sqrt(Lz*3),Zz,R2z)
    fprintf(' Dxy=%1.3fum^2/s, Lxy=%1.3fum, Axy=%1.2f, R2=%1.4f\n',Dxymicro,sqrt(Lxy*3/2),Zxymicro,R2xy)
    fprintf(' Dxyz=%1.3fum^2/s, Lxyz=%1.3fum, Ar=%1.2f, R2=%1.4f\n',Dr,sqrt(Lr),Zr,R2r)
end

%% VACF
XX=Xwalk;
YY=Ywalk;
ZZ=Zwalk;
ranmsdpsec=(1:sz(2))/frate;
if printtxt~=0
%     disp('VACF')
end
ranmsdf=50;
vacfdifco

%% PLOTS
% Plot MSD/VACF xyz
if plotMSDVACF==1
    figure
    %     subplot(1,2,1)
    plot(ranmsdosec',[Xmsd',Ymsd',Zmsd',XYmsd',XYZmsd']); hold on
    plot(fitX,':r');
    plot(fitY,':r')
    plot(fitZ,':r')
    plot(fitXY,':r')
    plot(fitR,':r')
    grid on
    ylabel(['MSD (',char(181),'m^2)'])
    xlabel('Lag time (s)')
    xlim([0 max(ranmsdosec)])
    axis square
    %     figure
    %     subplot(1,2,2)
    %     Tx=T1(1:FA1);
    %     [VACFplotxy,V1,S1]=plotyy(Tx,[VACFx(1:FA1),VACFy(1:FA1),VACFz(1:FA1),VACFxyz(1:FA1)],Tx,[DvacfxvsT2(1:FA1),DvacfyvsT2(1:FA1),DvacfzvsT2(1:FA1),DvacfxyzvsT2(1:FA1)]);
    %     set(get(VACFplotxy(1),'ylabel'),'string','VACF')
    %     %set(VACFplot(1),'ylim',[min(VACF) max(VACF)])
    %     set(get(VACFplotxy(2),'ylabel'),'string',['D (',char(181),'m^2/s)'])
    %     %set(VACFplot(2),'ylim',[0 2.5])
    %     xlabel('Time (s)')
    %     set(VACFplotxy(1),'xlim',[0 T1(FA1)])
    %     set(VACFplotxy(2),'xlim',[0 T1(FA1)])
    %     axis(VACFplotxy,'square')
    %     linkaxes([VACFplotxy(1),VACFplotxy(2)],'x')
end

XX=Xwalk;
YY=Ywalk;
ZZ=Zwalk;
if plot2D3D==1
    figure
    if reflectz==1
        plot3(XX,YY,[ZZ ZwalkNorm]);
    else
        plot3(XX,YY,ZZ);
    end
    xlabel(['x (',char(181),'m)'])
    ylabel(['y (',char(181),'m)'])
    zlabel(['z (',char(181),'m)'])
    grid on
    axis square
    set(gca,'dataaspectratio',[1 1 1]) % square pixels
end

if saveXYZwalkTs==1
    save('rwalk-trk.mat','Xwalk','Ywalk','Zwalk','Ts')
end

% aaa=sqrt((XX(2:end)-XX(1:end-1)).^2+(YY(2:end)-YY(1:end-1)).^2+(ZZ(2:end)-ZZ(1:end-1)).^2);

% figure
% [VACFplotxy,V1,S1]=plotyy(Tx,VACFxyz(1:ranmsdf+1),Tx,[DvacfxvsT2(1:ranmsdf+1),DvacfyvsT2(1:ranmsdf+1),DvacfzvsT2(1:ranmsdf+1),DvacfxyzvsT2(1:ranmsdf+1)]);
% set(get(VACFplotxy(1),'ylabel'),'string','VACF')
% %set(VACFplot(1),'ylim',[min(VACF) max(VACF)])
% set(get(VACFplotxy(2),'ylabel'),'string',['D (',char(181),'m^2/s)'])
% %set(VACFplot(2),'ylim',[0 2.5])
% xlabel('Time (s)')
% set(VACFplotxy(1),'xlim',[0 T1(ranmsdf+1)])
% set(VACFplotxy(2),'xlim',[0 T1(ranmsdf+1)])
% axis(VACFplotxy,'square')
% linkaxes([VACFplotxy(1),VACFplotxy(2)],'x')