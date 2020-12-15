clear all;
close all;

addpath('../tools');
load ../data/h_eig;
Platter = A;
Spin    = E;

%% General Setup

Bands = logspace(log10(100),log10(20000),25);
Fs    = 48000;


%% Calculate the RAPID Corrections
%  First load the data for the RAPID microphone waving

Files = {  '../data/RAPID_3.wav',
           '../data/RAPID_5.wav',
           '../data/RAPID_6.wav',
           '../data/RAPID_7.wav',
           '../data/RAPID_8.wav' };

for (f=1:length(Files))
    fprintf('%s\n',Files{f});
    [x fs] = audioread(Files{f}); 
    if (Fs ~= fs) x = resample(x,Fs,fs); end;
    H(:,:,f) = Response(x,Fs,[],Bands,[]);
end;

H_rapid = 10*log10(H+1E-4);

%% Figure showing the median of the microphone cluster for the three tries.
%  This is sort of a sense of how consistent the field and movement was in the
%  the sense of the three takes being repeatable.  Note that the median is
%  further subtracted in the calculation of rapid, which further reduces
%  bias
figure('name','Median Consistency','position',[100 100 640 480]);
axes('position',[.1 .15 .88 .83]);
semilogx(0,0,'b-','linewidth',2); hold on; set(gca,'fontsize',12);
semilogx(0,0,'m--','linewidth',2);
semilogx(Bands,squeeze(median(H_rapid,2)-mean(H_rapid(:))),'b-','linewidth',2); 
semilogx(Bands,squeeze(median(H_rapid,2)-mean(median(H_rapid,2),3)),'m-.','linewidth',2);
grid on;
axis([100 20000 -15 15]);
%axis([100 20000 -1 1]);
set(gca,'YTick',[-15 -10 -5 -1 0 1 5 10 15]);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude (dB)');
legend('Measured Specta','Spectra Difference to Mean');
print -dpng Fig04_RAPID_Median_Repeatability.png


%% Figure Plucking out a couple of individual microphones against the median
%  of the full set of microphones - this is our 'correction' according to
%  the RAPID approach.
figure('name','RAPID Delta to Median','position',[100 100 640 480]);
axes('position',[.1 .15 .88 .83]);
semilogx(0,0,'b-','Linewidth',2); hold on;
plot(0,0,'m--','Linewidth',2);
plot(0,0,'r-.','Linewidth',2);
semilogx(Bands,squeeze(H_rapid(:,22,:)-median(H_rapid,2)),'b-','linewidth',2); 
semilogx(Bands,squeeze(H_rapid(:,28,:)-median(H_rapid,2)),'m-.','linewidth',2); 
semilogx(Bands,squeeze(H_rapid(:,31,:)-median(H_rapid,2)),'r--','linewidth',2); 
grid on;
axis([100 20000 -6 2]);
xlabel('Frequency (Hz)');
ylabel('Mic Offset (dB)');
legend('Mic 22','Mic 28','Mic 31');
print -dpng Fig11_RAPID_Mic_Repeatability.png


%% Now plucking out a couple of responses against the median for the polar turntable plot
%  Note here that we add the energy of the responses weighted by the angle
%  density, which depends on the spin axis relative to the direction of the
%  speaker, or essentially the platter angle.
figure('name','Turntable Delta to Median','position',[100 100 640 480]);
H_polar = Response(h_eig,Fs,.2,Bands,[]);
H_polar = sum(sum(H_polar(:,:,1:end-1,:) .* permute(sin(Platter/180*pi),[1 3 4 2]),3),4);
H_polar = 10*log10(H_polar);
axes('position',[.1 .15 .88 .83]);
semilogx(0,0,'b-','Linewidth',2); hold on;
plot(0,0,'m--','Linewidth',2);
plot(0,0,'r-.','Linewidth',2);
semilogx(Bands,H_polar(:,22)-median(H_polar,2),'b-','linewidth',2); 
semilogx(Bands,H_polar(:,28)-median(H_polar,2),'m-.','linewidth',2); 
semilogx(Bands,H_polar(:,31)-median(H_polar,2),'r--','linewidth',2); 
grid on;
axis([100 20000 -6 2]);
xlabel('Frequency (Hz)');
ylabel('Mic Offset (dB)');
legend('Mic 22','Mic 28','Mic 31');
print -dpng Fig12_Polar_Mic_Diffuse.png


%% Take a large set of mics and plot the average of the hand waving
%  against the ground truth from the turntable.
%
%  Note that arguably, the turntable has more systematic bias than the
%  RAPID due to the scatter of the turn table.
%
figure('name','Turntable vs Delta','position',[100 100 640 480]);
M = [ 3 5 8 10 13 14 17 22 27 29 31 32 ];
for (m=1:length(M))
    r = 3-floor((m-1)/3);
    c = mod(m-1,3);
    axes('position',[.06+.31*c .06+.23*r .30 .22]);
    % calibrate against the median 
    rapdidMed = squeeze(H_rapid(:,M(m),:)-median(H_rapid(:,:,:),2));
    polarMed = H_polar(:,M(m))-median(H_polar,2);
    
    %interpolate for pretty plots
    rapidInterp = interp1([0 25 Bands 3e4],[[1; 1;].*rapdidMed(1,:);  rapdidMed(:,:); rapdidMed(end,:)],(0:Fs/2),'pchip');
    polarInterp = interp1([0 25 Bands 3e4],[[1; 1;].*polarMed(1);  polarMed(:,:); polarMed(end,:)],(0:Fs/2),'pchip');
    
    h1 = semilogx((0:Fs/2),rapidInterp,'b:','LineWidth',2); hold on;
    h2 = semilogx((0:Fs/2),polarInterp,'r-','LineWidth',2.5); hold on;
    grid on;
    axis([100 20000 -1.5 1.5]);
    set(gca,'fontsize',10,'TickLabelInterpreter','latex');
    if (r>0) set(gca,'xticklabel',[]); end;
    if (r==0) set(gca,'xtick',[1000 10000]); set(gca,'xticklabel',{'$1$ kHz', '$10$ kHz'}); end;
    if (c>0) set(gca,'yticklabel',[]); end;
    if (c==0) set(gca,'yticklabel',[-1 0 1]); set(gca,'yticklabel',{'$-1$ dB','$0$ dB','$1$ dB'}); end;
    text(120,-1.2,['$\textrm{Mic}$ ' num2str(M(m))],'interpreter','latex'); set(gca, 'FontName', 'Times New Roman')
end;
%
legend1 = legend([h1(end) h2], 'RAPID', 'Turn Table');
set(legend1,'Position',[0.52 0.99 0 0],'Orientation','horizontal','FontSize',12);  
legend boxoff  
print -dpng Fig05_Polar_vsRAPIDs.png

