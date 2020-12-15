%% Create an array of beam plots
function BeamPlots(X, F, A, Lim, Type, C, CTick, FL, CL, FT)

if (nargin<5) Lim = [min(X(:)) max(X(:))]; end;
if (nargin<5) Type = 1; end;
if (nargin<6) C = parula;  end;
if (nargin<7) CTick = Lim(1)+(Lim(2)-Lim(1))*[0 .22 .4 .6 .8 1];  end;
if (nargin<8) FL = 'Hz'; end;
if (nargin<9) CL = 'dB'; end;
if (nargin<10) FT = 7; end;

N = size(X,3);
Rows = round(sqrt(N));
Cols = ceil(N/Rows);

clf;
colormap(C);
axes('position',[ .04 .06 .92 .02 ]);   % Add a color legend
image(repmat(1:length(C),5,1));
set(gca,'fontsize',10);
set(gca,'ytick',[]);
set(gca,'xtick',(length(C)-1)*(CTick-Lim(1))/(Lim(2)-Lim(1))+.5);
set(gca,'xticklabel',arrayfun(@(x)sprintf('%.0fdB',x),CTick,'UniformOutput',false));


w = (.92-(Cols-1)*.02)/Cols;  h = (.88-(Rows-1)*.02)/Rows;  g = .02;

if (Type==2) [A I] = sort(A); X = X(:,I,:); end;

for (n=1:N)
     r = Rows - 1 - floor((n-1)/Cols);
     c = mod(n-1,Cols);
     axes('position',[ .04+c*(w+g) .1+r*(h+g) w h ]);
     switch (Type)
        case 1
            RingPlot(X(:,:,n)',A,.25,1,2); hold on;
            caxis(Lim);
            set(gca,'PlotBoxAspectRatio',[1 1 1]);
            set(gca,'xtick',[]); set(gca,'ytick',[]);
            set(gca,'FontSize',16); 
            set(gca,'XColor', 'none','YColor','none');
            set(gca,'color','none');
            plot([.25 1],[0 0],'k-','linewidth',1); 
            text(0,0,sprintf('%d',n),'horizontalalignment','center','verticalalignment','middle');
            if (n==1)
                text(0, .20,sprintf('%.0f %s',F(1),FL)  ,'HorizontalAlignment','center','VerticalAlignment','middle');
                text(0,1.0,sprintf('%.0f %s',F(end),FL),'HorizontalAlignment','center','VerticalAlignment','middle');
                text(0.5,.05,'0 deg');
            end;
        case 2
            surf(A,(1:length(F))-1,X(:,:,n),'CDataMapping','scaled'); view(2); shading interp; hold on;
            plot3([0 0],[1 length(F)],max(X(:))*[1 1],'k');
            caxis([Lim(1) Lim(2)]);
            set(gca,'xtick',[]); set(gca,'ytick',[]);
            xlim([A(1) A(end)]); ylim([0 length(F)-1]);
            set(gca,'XColor', 'none','YColor','none');
            set(gca,'color','none');
            if (n==1)
                text(0,0.01*length(F),max(X(:)),sprintf('0 deg'),'HorizontalAlignment','center','color',.7*[1 1 1]);
                text(max(A),0.01*length(F),max(X(:)),sprintf('%.0f deg ',max(A)),'HorizontalAlignment','right','color',.7*[1 1 1]);
                text(min(A),0.99*length(F)-1,max(X(:)),sprintf('  %.0f %s',F(end),FL),'HorizontalAlignment','left','verticalalignment','top','color',.7*[1 1 1]);
                for(f = (length(F)-1)/(FT-1)+1:(length(F)-1)/(FT-1):length(F)-1)
                    text(min(A),f-.5,max(X(:)),sprintf('  %.0f %s',F(ceil(f)),FL),'HorizontalAlignment','left','color',.7*[1 1 1]);
                end;
                text(min(A),0.01*length(F),max(X(:)),sprintf('  %.0f %s',F(1),FL),'HorizontalAlignment','left','color',.7*[1 1 1]);
            end;
    end;
end;
        