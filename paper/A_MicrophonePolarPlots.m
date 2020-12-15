clear all;

FigSize = [ 100 100 640 480 ];
addpath('../tools');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DBX MICROPHONE CALIBRATION DATA
%  
%  Load in the detais for the DBX micorphone that we will be using.
%  This includes the calibration data of this mic against 5 other reference
%  mics for the on axis response, and then 2 axis turntable data for off 
%  corrections.  This is important as we will not have the microphone pointing
%  on axis at all of the speakers, and we want to correct for the point response.
%

if (exist('../data/DBX_OffAxis.mat'))
    load ../data/DBX_OffAxis.mat;
else
    Bands = logspace(log10(100),log10(20000),25);

    load('../data/DBX.cal');
    DBX     = spline(DBX(:,1),DBX(:,2),Bands)';

    E = [ 0:10:180 ]';
    A = [ 0:10:360 ]';
    Delay  = 0.085;
    Gain   = 0.1;
    Length = 0.1;
    SWEEP_20200609_48000_2000
    X = zeros(Length*Fs,length(E),length(A));
    for (e=1:length(E))
        for (a=1:length(A))
            file = sprintf('../data/DBX/20200609_COMOTTABLE_DBX_E%03.0f_A%03.0f.wav',E(e),A(a));
            fprintf('Loading %s\n',file);
            x = audioread(file);
            x = (1/Gain)*Convolve(x,SweepI);
            x = x(Offsets(1)+Delay*Fs+(1:size(X,1)),:);
            X(:,e,a) = x;
        end;
    end;

    DBX_Hpow = Response(X,Fs,Length,Bands,[0.001 Length/2]);
    DBX_H = 10*log10(DBX_Hpow+.0000001);
    DBX_E = mean(DBX_H(:,:,:),3) - mean(DBX_H(:,1,:),3);        % Average around azimuth, correct on axis
    DBX_E(1:4,:) = 0;                                   % Remove LF estimation noise
    DBX_E = DBX_E + DBX;
    DBX_H = DBX_H;
    DBX_X = X;
    save ../data/DBX_OffAxis Bands E DBX_E DBX DBX_H DBX_Hpow DBX_X Fs
end;

%% Figure Reference microphone with angle

% Calculate the polar eq
DBX_Hdiff = 10*log10(sum(sum(DBX_Hpow(:,:,1:end-1) .* permute(sin(E/180*pi),[2 1 3]),2),3)+.0000001);
DBX_Hdiff = DBX_Hdiff - 10*log10(sum(sum(repmat(DBX_Hpow(:,1,1:end-1),1,length(E),1) .* permute(sin(E/180*pi),[2 1 3]),2),3)+.0000001);

Bands = logspace(log10(100),log10(20000),25);


figure('name','DBX Off Axis','position',FigSize);
axes('position',[.1 .75 .88 .23]); set(gca,'fontsize',12);
semilogx(Bands,DBX_Hdiff,'b-','linewidth',2); hold on;
plot(0,0,'color',[1 1 1]);
plot(0,0,'color',[1 1 1]);
axis([100 20000 -5 1]);
set(gca,'YTick',[-15 -10 -5 -1 0 1]);
set(gca,'xtickLabel',[]);
grid on;
legend('Diffuse Response','Weighted sum of','Power over Elevation',...
    'location','southwest','Interpreter','latex','fontsize',12);

axes('position',[.1 .15 .88 .57]); set(gca,'fontsize',12);
semilogx(0,0,'b-','linewidth',2); hold on;
semilogx(0,0,'g:','linewidth',2);
semilogx(0,0,'r--','linewidth',2);
semilogx(0,0,'m-','linewidth',2);
semilogx(0,0,'c:','linewidth',2);
semilogx(0,0,'y--','linewidth',2);
semilogx(0,0,'k:','linewidth',2);
semilogx(Bands,squeeze(DBX_H(:,1,1:3:end)-mean(DBX_H(:,1,:),3)),'b-'); 
semilogx(Bands,squeeze(DBX_H(:,4,1:3:end)-mean(DBX_H(:,1,:),3)),'g:');
semilogx(Bands,squeeze(DBX_H(:,7,1:3:end)-mean(DBX_H(:,1,:),3)),'r--'); 
semilogx(Bands,squeeze(DBX_H(:,10,1:3:end)-mean(DBX_H(:,1,:),3)),'m-'); 
semilogx(Bands,squeeze(DBX_H(:,13,1:3:end)-mean(DBX_H(:,1,:),3)),'c:'); 
semilogx(Bands,squeeze(DBX_H(:,16,1:3:end)-mean(DBX_H(:,1,:),3)),'y--'); 
semilogx(Bands,squeeze(DBX_H(:,19,1:3:end)-mean(DBX_H(:,1,:),3)),'k:'); 
legend('$0^{\circ}$','$30^{\circ}$','$60^{\circ}$','$90^{\circ}$','$120^{\circ}$',...
    '$150^{\circ}$','$180^{\circ}$','location','southwest','Interpreter','latex','fontsize',12);
grid on;
axis([100 20000 -15 1]);
set(gca,'YTick',[-15 -10 -5 -1 0 1]);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude (dB)');
set(gca,'fontsize',10,'TickLabelInterpreter','latex');
set(gca, 'FontName', 'Times New Roman')
print -dpng Fig03_POLAR_DBX_Diffuse.png

%% Do a polar plot and contour plot of the reference mic

Bands = logspace(log10(2000),log10(16000),13);

Hpow = Response(DBX_X,Fs,0.1,Bands,[0.001 0.1/2]);
H = 10*log10(Hpow+.0000001);

E = [ 0:10:180 ]';
G = H(:,:,:) - mean(H(:,1,:),3);
G = [ G(:,end:-1:2,19) G(:,:,1) ];
Eb = [-flipud(E(2:end)); E];

figure('name','DBX Flame','position',FigSize);
axes('position',[.1 .75 .88 .23]); set(gca,'fontsize',12);
BeamPlots(G,Bands,Eb,[-5 1],2,jet,[-5 -4 -3 -2 -1 0 1],'Hz','dB',4);
hold on;
c = colormap; 
cm3 = 0.8*c(round(length(colormap)*17/30),:);
contour(Eb,0:length(Bands)-1,G,[-5 -4 -3 -2 -1],'linecolor','k');

print -dpng Fig01_POLAR_DBX_Flame.png

%% Do the same plot, but for one mic of the Eigen Mic
load ../data/h_eig;        % Turntable response of EigenMic already deconvolved
Platter = A;
Spin    = E;

Mic     = 2;
P       = 6;

H_polar = 10*log10(Response(h_eig,Fs,.2,Bands,[])+.00001);
G = squeeze(H_polar(:,Mic,:,P));
G = G - G(:,3);

Spin = -180:15:180;
G = G(:,[18:25 1:17]);

figure('name','EIG Flame','position',FigSize);
axes('position',[.1 .75 .88 .23]); set(gca,'fontsize',12);
BeamPlots(G,Bands,Spin,[-5 1],2,jet,[-5 -4 -3 -2 -1 0 1],'Hz','dB',4);
hold on;
c = colormap; 
cm3 = 0.8*c(round(length(colormap)*17/30),:);
contour(Spin,0:length(Bands)-1,G,[-5 -4 -3 -2 -1],'linecolor','k');
print -dpng Fig02_POLAR_EIG_Flame.png
