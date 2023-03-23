% Velocity autocorrelation function by straigh forward and fft methods for
% calculating diffusion coefficient
% Teeranan, 2015

plotvacf=0;
pxlsize = 65.85; %65.85nm/pixel
Fst=5; % Start frame number to integrate

FA1=FA+1;

tsec=ranmsdpsec';

% Vx=diff(XX)/mean(diff(tsec));
% Vy=diff(YY)/mean(diff(tsec));
% % Vz=(1/(2*2e-5)*diff(WW*1000))/1000/mean(diff(tsec));
% Vz=diff(ZZ)/mean(diff(tsec));
% 
% Nx=length(Vx);
% 
% Vx2=[0;Vx];
% Vy2=[0;Vy];
% Vz2=[0;Vz];
% Vxyz=[Vx2,Vy2,Vz2];
% VACFxyzT=zeros(3,Nx)';
% VACFfftxy=zeros(2,Nx)';

% tic
% Vxyz=[0 0 0;diff(XYZwalkt)]*frate;
% VACFxyzT=zeros(FA,3);
% for i=1:3
%     V=Vxyz(:,i);
%     Nv=length(V);
%     
% %     (1) Straightforward code for autocorrelation
%     for j=1:FA
%         for k=1:Nv-j+1
%             VACFxyzT(j,i)=VACFxyzT(j,i)+V(k)*V(j+k-1);
%         end
%     end
%     VACFxyzT(:,i)=VACFxyzT(:,i)/Nv;
%     
% %     % (2) FFT for autocorrelation
% %     nfft=2^nextpow2(2*Nv-1);
% %     Cfft=real(ifft(fft(V,nfft).*conj(fft(V,nfft))));
% %     VACFxyzT(1:FA,i)=Cfft(1:FA)/Nv;
% end
% vacftime(end+1)=toc;
tic
Vxyz=[0 0 0;diff(XYZwalkt)]*frate;
VACFxyzT=zeros(FA1,3);
Nv=length(Vxyz);
norm=(Nv:-1:Nv-FA)'; % Normalizer
for i=1:3
    V=Vxyz(:,i);
% %     (1) Straightforward code for autocorrelation
%     for j=1:FA
%         for k=1:Nv-j+1
%             VACFxyzT(j,i)=VACFxyzT(j,i)+V(k)*V(j+k-1);
%         end
%     end
%     VACFxyzT(:,i)=VACFxyzT(:,i)./norm;
    
    % (2) FFT for autocorrelation
    nfft=2^nextpow2(2*Nv-1);
    Cfft=real(ifft(fft(V,nfft).*conj(fft(V,nfft))));
    VACFxyzT(1:FA1,i)=Cfft(1:FA1)./norm;
end

if exist('vacfffttime','var')
    vacfffttime(end+1)=toc;
end
VACFxy=VACFxyzT(:,1)+VACFxyzT(:,2);
VACFx=VACFxyzT(:,1);
VACFy=VACFxyzT(:,2);
VACFz=VACFxyzT(:,3);
VACFxyz=VACFxy+VACFz;

ts=[0;tsec(1:Nv-1)];
T1=ts(1:length(VACFxy));

sumareaxy=[];
for i=2:FA1
    sumareaxy(i,1)=trapz(T1(1:i),VACFxy(1:i));
end
DvacfxyvsT2=sumareaxy(1:end)/2;

sumareax=[];
for i=2:FA1
    sumareax(i,1)=trapz(T1(1:i),VACFx(1:i));
end
DvacfxvsT2=sumareax(1:end);

sumareay=[];
for i=2:FA1
    sumareay(i,1)=trapz(T1(1:i),VACFy(1:i));
end
DvacfyvsT2=sumareay(1:end);

sumareaz=[];
for i=2:FA1
    sumareaz(i,1)=trapz(T1(1:i),VACFz(1:i));
end
DvacfzvsT2=sumareaz(1:end);

sumareaxyz=[];
for i=2:FA1
    sumareaxyz(i,1)=trapz(T1(1:i),VACFxyz(1:i));
end
DvacfxyzvsT2=sumareaxyz(1:end)/3;

T2=ts(2:end);

if plotvacf==1
    Tx=T1(1:FA1);
    figure
    subplot(1,2,1)
    [VACFplotxy,V1,S1]=plotyy(Tx,[VACFx(1:FA1),VACFy(1:FA1),VACFxy(1:FA1)],Tx,[DvacfxvsT2(1:FA1),DvacfyvsT2(1:FA1),DvacfxyvsT2(1:FA1)]);
    set(get(VACFplotxy(1),'ylabel'),'string','VACFxy')
    %set(VACFplot(1),'ylim',[min(VACF) max(VACF)])
    set(get(VACFplotxy(2),'ylabel'),'string',['D (',char(181),'m^2/s)'])
    %set(VACFplot(2),'ylim',[0 2.5])
    xlabel('Time (s)')
    set(VACFplotxy(1),'xlim',[0 T1(FA1)])
    set(VACFplotxy(2),'xlim',[0 T1(FA1)])
    axis(VACFplotxy,'square')
    linkaxes([VACFplotxy(1),VACFplotxy(2)],'x')

    subplot(1,2,2)
    [VACFplotz,V2,S2]=plotyy(Tx,VACFz(1:FA1),Tx,DvacfzvsT2(1:FA1));
    set(get(VACFplotz(1),'ylabel'),'string','VACFz')
    %set(VACFplot(1),'ylim',[min(VACF) max(VACF)])
    set(get(VACFplotz(2),'ylabel'),'string',['D (',char(181),'m^2/s)'])
    %set(VACFplot(2),'ylim',[0 2.5])
    xlabel('Time (s)')    
    set(VACFplotz(1),'xlim',[0 T1(FA1)])
    set(VACFplotz(2),'xlim',[0 T1(FA1)])
    axis(VACFplotz,'square')
    linkaxes([VACFplotz(1),VACFplotz(2)],'x')

end

% DvacfxyAll=mean(DvacfxyvsT2(:));
% fprintf(' Dxy(5toEnd)=%1.3fum^2/s\n',DvacfxyAll)
Dvacfx=mean(DvacfxvsT2(Fst:FA1));
fprintf(' Dx=%1.3fum^2/s\n',Dvacfx)
Dvacfy=mean(DvacfyvsT2(Fst:FA1));
fprintf(' Dy=%1.3fum^2/s\n',Dvacfy)
Dvacfz=mean(DvacfzvsT2(Fst:FA1));
fprintf(' Dz=%1.3fum^2/s\n',Dvacfz)
% DvacfzAll=mean(DvacfzvsT2(:));
% fprintf(' Dz(5toEnd)=%1.3fum^2/s\n',DvacfzAll)
Dvacfxy=mean(DvacfxyvsT2(Fst:FA1));
fprintf(' Dxy=%1.3fum^2/s\n',Dvacfxy)
Dvacfxyz=mean(DvacfxyzvsT2(Fst:FA1));
fprintf(' Dxyz=%1.3fum^2/s\n',Dvacfxyz)

% Mass from VACF(0)
% 1.38e-23*296/(VACFxyzT(1,2)*Nv*1e-12*1e-6)