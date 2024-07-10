function Plot_Bump_Examples(DarkData, ttt, Seg)
% Plot 40s examples of the bump for specified trials and segments

PlotWindow=40; % seconds

Time = DarkData.TwoP_Times(ttt,:);
Time = Time - Time(1);

% Get indices for this segment
starttemp = (Seg-1)*PlotWindow;
stoptemp  = starttemp + PlotWindow;
logica    = Time >= starttemp & Time <= stoptemp;

% Get data to plot
PlotData      =  DarkData.ROI_Data_Dfof(:,logica,ttt);
WedgeCenters  =  DarkData.WedgeCenters;
TimeTemp      =  Time(logica);
FlyHD         =  DarkData.FlyPosition_w2p_aligned(ttt, logica);
BumpPos       =  DarkData.BumpPosition(ttt, logica);
BumpR         =  DarkData.BumpR(ttt, logica);

% Compute the Fly HD from the Fly AV
if ttt<13
    FlyAV    =  DarkData.FlyRotate_w2p(ttt,:);
    dT       =  1/DarkData.TwoP_Fs(1);
    FlyHD_2  =  cumsum(-FlyAV*dT)/360*2*pi;
    FlyHD_2  =  mod(FlyHD_2, 2*pi)-pi;
    FlyHD    =  FlyHD_2(logica);
end

% Calculate offset
R_Thresh=prctile(BumpR,80);
FlyPositionsStrong=FlyHD(logical(BumpPos(BumpR>R_Thresh)));
BumpPositionsStrong=BumpPos(logical(BumpPos(BumpR>R_Thresh)));
AngleDiff=FlyPositionsStrong-BumpPositionsStrong;
AngleDiff(AngleDiff>=pi)=AngleDiff(AngleDiff>=pi)-2*pi;
AngleDiff(AngleDiff<=-pi)=AngleDiff(AngleDiff<=-pi)+2*pi;
Offset=circ_mean(AngleDiff');

% Compute aligned arena position
FlyHD=FlyHD-Offset;
FlyHD(FlyHD<-pi)=FlyHD(FlyHD<-pi)+2*pi;
FlyHD(FlyHD>pi)=FlyHD(FlyHD>pi)-2*pi;

% Get rid of wraparound artifacts in FlyHD
FlyHD_ToPlot = FlyHD;
wrapped = abs(diff(FlyHD)) > 1.5*pi;
FlyHD_ToPlot(wrapped) = nan;

% Plot the bump data
figure;
hold on;
imagesc(TimeTemp, WedgeCenters,PlotData);
map = brewermap(256,'Blues');
colormap(map)
colormax=ceil(prctile(PlotData(:),97));
cb=colorbar;
clim=get(gca,'clim');
set(gca,'clim',[0 floor(clim(2)/10)*10])
set(cb,'Location','eastoutside')
set(gca,'clim',[0 colormax])
xlim([TimeTemp(1) TimeTemp(end)])
ylim([-pi pi])
set(gca,'ytick',-pi:pi/4:pi)
set(gca,'xtick', floor(TimeTemp(1):10:TimeTemp(end)))
ylabel('Wedge Pos. (rad)')
plot(TimeTemp, FlyHD_ToPlot, 'r');
xlabel('Time (s)');
axis equal;
axis tight;
end

