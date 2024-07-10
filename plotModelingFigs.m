function plotModelingFigs(plotInds)

%---------------------------- KEY FOR PLOTINDS ----------------------------%
%1: plotConnMat_Schematic                   % Fig 2A    
%2: plotOptJE                               % Fig 2D
%3: plotEnergyHeatmaps                      % Fig 2E
%4: plotBumpTraj_driven                     % Fig 2F    
%5: plotBumpTraj_noinput                    % Fig 2F-G  
%6: plotConnMat_VarJ                        % Fig 2H

% FIGURE 3
%7:  plotBumpTraj_noinput_regimes           % Fig 3C
%8:  plotRegimeRates                        % Fig 3D;
%9:  plotSimplifiedEnergy_noinput           % Fig 3E-F
%10: plotBumpTraj_noinput_regimes_varJE     % Fig 3G; ED Fig 8G
%11: plotNetDriftSpeed                      % Fig 3H
%12: plotSimplifiedEnergy_driven            % Fig 3I-J
%13: plotBumpTraj_driven_regimes_varJE      % Fig 3K; ED Fig 8G
%14: plotThreshVel                          % Fig 3L
%15: plotLinearity                          % Fig 3L

% FIGURE 4
%16: plotParameterRobustness                % Fig 4A-C; ED Fig 9A 
%17: plotNoiseRobustness                    % Fig 4D-F, ED Fig 9B

% ADDITIONAL ED FIGURE PANELS
%18: plotPhaseDiagram                       % ED Fig 3A
%19: plotFixedPointConditions               % ED Fig 3B,C
%20: plotLinearizationEvals                 % ED Fig 3D-I
%21: plotHessianEvals                       % ED Fig 4
%22: plotLeadingEvals_noinput               % ED Fig 6
%23: plotLeadingEvals_smallinput            % ED Fig 7
%24: plotRegimeSpeeds                       % ED Fig 8B
%25: plotFixedPoints_regimes_varV           % ED Fig 8C,E                                                 
%26: plotFixedPoints_regimes_belowThreshV   % ED Fig 8D
%27: plotFixedPoints_regimes_aboveThreshV   % ED Fig 8F

%-------------------------------------------------------------------------%


%------------------------ plotting parameters ----------------------------%
cT = [0,149,163]./255;      % stable regime
cO = [255,135,0]./255;      % unstable regime
cW = [1,1,1];               % white 
cK = [0,0,0];               % black
cR = [255,0,54]./255;       % red color for different values of JE
cmapN = [[72,13,46];        % colormap for plotting N-dependence
    [107,16, 56 ];
    [130,40, 83 ];
    [169,91, 138];
    [206,146,190];
    [240,199,241];
    [245,220,245];
    [253,245,253]]./255;
cmapN = flipud(cmapN);

lw = 2;                     % linewidth
ms = 12;                    % marker size
fs = 16;                    % font size
alphaR = 0.5;               % opacity

%------------------------ simulations parameters -------------------------%
N  = 6;                                                 % number of computational units
C0 = 1;                                                 % constant feedforward input
A  = 0.2;                                               % approx amplitude to aim for (au, used in computation of JI)
tau = 0.1;                                              % neural time constant
del = pi/2;                                             % phase shift for velocity-dependent weights

th = linspace(0, 2*pi, N+1);  th(end) = [];             % preferred feature space
psiSamp = linspace(0,2*pi,101); psiSamp(end) = [];      % sampled bump orientations
dth     = th(2)-th(1);                                  % angular separation between orientations

[Ws,Wa,~] = getConnectivity(N,3,A,C0,psiSamp,del,0);    % symmetric & asymmetric components of connectivity matrices
W = getFullConnectivity(Ws,Wa,0,tau);                   % full connectivity matrix in absence of velocity input

[JEopt,~,NactOpt] = computeOptJE(N);                    % optimal excitation for N = 6

%----------------------------- for plotting ------------------------------%

% time increments to use for plotting
dt    = .01;                                    % time increment
tmax1 = 3;                                      % short max time interval                  
tmax2 = 6;                                      % long max time interval
tV = 0:dt:tmax2;                                % velocity scaling

% range of JE to use for plotting (1/JEmax < 1/JE < 1/JEmin)
JEmin = JEopt(1);
JEmax = JEopt(2);

% suboptimal values of JE
JEtmp   = 1./linspace(1./JEmax,1./JEmin,25);
JEsampN = JEtmp([2,5,13,21,24]);                % JE = [3.89, 3.6, 3, 2.57, 2.44]
JEsamp0 = JEsampN(3);                           % JE = 3 (used for most comparisons)     
JEsampE = JEtmp([2,5,9,13,17,21,24]);           % used to illustrate energy landscapes


%values to use for comparison with analytics
JEsim = 1./linspace(1./JEmax,1./JEmin,25);
JEsim = JEsim(2:end-1);

%specific values of JE to use across all JE space ('W' = wide sampling)
JEsampW = JEopt(1);
OsampW  = 1;   
NactW   = NactOpt(1);
for i=1:numel(JEopt)-1
    JEtmp   = fliplr(1./linspace(1./JEopt(i+1),1./JEopt(i),3));
    JEsampW = [JEsampW,JEtmp(2:3)];
    OsampW  = [OsampW,[0,1]];
    NactW   = [NactW,NactOpt(i:i+1)];
end
JEsampW = fliplr(JEsampW);
OsampW  = fliplr(OsampW);
NactW   = fliplr(NactW);


%initial values of psi to use for bump drift
ntraj = 6;
psi0tmp = linspace(0,pi/N,ntraj);
psi0tmp = psi0tmp(1:end-1);
psi0    = [];
for i=1:2*N
    psi0 = [psi0,psi0tmp+(i-1)*(pi/N)];
end

%perturbation for drift
dpsi = 0.01;

%for velocity-driven plots
[~,~,tauU,tauS] = computeTimeConstants(JEsamp0,JEmin,JEmax,tau);
[~,dpsiS] = computeDelNact(tauU,tauS,N);
vmin0 = computeVmin(tauS,dpsiS);

%velocities used by default for plotting
vdesired_default  = .2:.2:2;                

%velocities used to illustrate bump trajectories for different JE
mdesired_bumpTraj_varJE = .1:.1:1;                

%velocites used to illustrate linearity of integration
mdesired_linearity = .8:.2:1.6;

%velocites used to illustrate bump trajectories below/above threshold velocity
mdesired_bumpTraj = linspace(0,vmin0,4);
dv = diff(mdesired_bumpTraj(1:2));
mdesired_bumpTraj = [mdesired_bumpTraj,mdesired_bumpTraj(end)+dv,mdesired_bumpTraj(end)+2*dv];

%velocities used to plot fixed points below and above threshold velocity
mdesired_FPbelowThresh = linspace(0,vmin0,9);               
mdesired_FPbelowThresh(1) = [];
          
mdesired_FPaboveThresh = linspace(0,3*vmin0,25);             
mdesired_FPaboveThresh(1) = [];

%--------------------------- build sim structures ------------------------%

plotVars.cK = cK;
plotVars.cW = cW;
plotVars.cO = cO;
plotVars.cT = cT;
plotVars.cR = cR;
plotVars.lw = lw;
plotVars.fs = fs;
plotVars.ms = ms;
plotVars.alphaR = alphaR;
plotVars.cmapN  = cmapN;

simVars.A   = A;
simVars.C0  = C0;
simVars.tau = tau;
simVars.dt  = dt;
simVars.del = del;
simVars.th  = th;
simVars.dth = dth;
simVars.ntraj   = ntraj;
simVars.dpsi    = dpsi;
simVars.psiSamp = psiSamp;

%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
if ismember(1,plotInds)
    % plot symmetric and asymmetric connectivity matrices for JEsamp0
    
    figure;set(gcf,'Position',[200 200 800 300],'color','w')
    
    [Ws,Wa] = getConnectivity(N,JEsamp0,A,C0,psiSamp,del);
    
    subplot(1,2,1);
    imagesc(Ws);colormap(gray); axis off
    title('symmetric connectivity')
    set(gca,'fontsize',fs);
    
    subplot(1,2,2);
    imagesc(Wa);colormap(gray); axis off
    title('asymmetric connectivity')
    set(gca,'fontsize',fs);
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(2,plotInds)
    % plot optimal JE vs N
    
    figure;set(gcf,'Position',[200 200 800 400],'color','w');hold on;
    
    %initialize
    Nall     = 4:22;
    JEinvTmp = nan(numel(Nall),max(Nall)-3);
    NtotTmp  = nan(numel(Nall),max(Nall)-3);
    NactTmp  = nan(numel(Nall),max(Nall)-3);
    NdiffTmp = nan(numel(Nall),max(Nall)-3);
    for i=Nall
        [~,JEinvTmp(i-3,1:i-3),NactTmp(i-3,1:i-3)] = computeOptJE(i);
        NtotTmp( i-3,1:i-3) = i;
        NdiffTmp(i-3,1:i-3) = NactTmp(i-3,1:i-3)-i/2; 
    end
    
    nd = unique(NdiffTmp(~isnan(NdiffTmp)));
    
    for i=1:numel(nd)
        inds = find(NdiffTmp == nd(i));
        if mod(2*nd(i),2)<.5
            plot(JEinvTmp(inds),NtotTmp(inds),'-ok' ,'markeredgecolor','none','markerfacecolor',cR,'markersize',ms);
        else
            plot(JEinvTmp(inds),NtotTmp(inds),'--sk','markeredgecolor','none','markerfacecolor',cR,'markersize',ms);
        end
    end

    ylim([4,20]);
    set(gca,'fontsize',fs,'TickDir','out')
    ylabel('number of neurons N')
    xlabel('excitation 1/J_E')
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(3,plotInds)
    % plot Energy Heatmaps
    
    figure;set(gcf,'Position',[200 200 2000 300],'color','w');
    
    %densely sample psi and thc
    nS = 500;
    thcD = linspace(pi/N, (N-1)*pi/N, nS+1);thcD = thcD';
    psiD = linspace(0,2*pi, nS+1);          psiD = psiD';
    
    %define color ranges
    yminOpt  = linspace(1.65,3.65,5);
    dcOpt    = fliplr([.0004,.0015,.0015,.005,.006]);

    
    nJ = numel(JEsampW);
    for i = 1:nJ
        
        [~,~,JI] = getConnectivity(N,JEsampW(i),A,C0,psiD,del);
        E = energyLandscape(JEsampW(i),JI,C0,N,psiD,thcD);
        
        %identify minima
        [~, inds] = min(E);
        minPsi = psiD;
        minThc = thcD(inds);
        
        %identify contours of constant frho (defines minima)
        frho = zeros(length(thcD),length(psiD));
        for p = 1:length(psiD)
            [f,~] = calculateFunctions([psiD(p)*ones(size(thcD)),thcD],N);
            frho(:,p) = f(:,2);
        end
        Mc = contourc(psiD, thcD, frho, [1/JEsampW(i) 1/JEsampW(i)]);
        psi = Mc(1,2:end)';
        thc = Mc(2,2:end)';
       

        %plot energy landscape
        subplot(1,nJ,i);
        imagesc(psiD, 2*thcD, E);colormap(gray);
        ylim( [yminOpt(i),  yminOpt(i)+pi/2 ]) 
        clim([min(min(E)), min(min(E))+dcOpt(i)])
        title(['1/J_E = 1/', num2str(JEsampW(i))])
        set(gca,'YDir','normal','fontsize',fs)
        axis off
        hold on
        
        if OsampW(i)<.5
            psiFP = [th 1/2*(th(1:end-1) + th(2:end)) 1/2*(th(end) + 2*pi) 2*pi];
            thcFP = interp1(psi,thc,psiFP);
            
            %plot energy minima
            plot(psiFP, 2*thcFP, 'wo','linewidth',2)
        else
            %plot minimum energy contour
            plot(psi, 2*thc, 'w','linewidth',2)
        end
        
    end
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(4,plotInds)
    % plot driven bump for different JE
    
    figure;set(gcf,'Position',[200 200 2000 800],'color','w');
    
    nJ = numel(JEsampW);
    
    psiInit = th(1).*ones(nJ);
    psiInit(OsampW<.5 & (1:nJ)<nJ/2) = (th(1)+th(2))./2;
   
    for j=1:nJ
        subplot(2,nJ,j);hold on;
        plotMultDrivenTrajAboveThresh(vdesired_default,JEsampW(j),N,tmax1,tmax2,simVars,plotVars,psiInit(j))
        
        subplot(2,nJ,j+nJ);hold on;
        plotMultDrivenTrajAboveThresh(vdesired_default,JEsampW(j),N,tmax1,tmax2,simVars,plotVars,psiInit(j))
        ylim([0,1])
        xlim([psiInit(j),psiInit(j)+pi/N])
        xticks([psiInit(j),psiInit(j)+pi/N])
        xticklabels({'0','\pi/N'}); 

    end  

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(5,plotInds)
    % plot drifting bump for different JE
    
    figure;set(gcf,'Position',[200 200 2000 800],'color','w');
     
    nJ   = numel(JEsampW);
    for j=1:nJ
        subplot(2,nJ,j);hold on;
        [~,~,bumpProfiles,indsbump] = plotMultDriftTraj(JEsampW(j),N,tmax1,tmax2,psi0,OsampW(j),0,simVars,plotVars);
        plot(indsbump,tmax1,'o')
        
        subplot(2,nJ,j+nJ);hold on;
        plotBumpProfiles(bumpProfiles,OsampW(j),simVars,plotVars);
         
    end
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(6,plotInds)
    % plot symmetric and asymmetric connectivity matrices for JEsampW
    
    figure;set(gcf,'Position',[200 200 1800 600],'color','w')
    
    nJ = numel(JEsampW);
    for i=1:nJ
        [Ws,~,~] = getConnectivity(N,JEsampW(i),A,C0,psiSamp,del);
        dW = (max(Ws(:))-min(Ws(:)))/2;
        subplot(2,nJ,i);
        imagesc(Ws-min(Ws(:))-dW);colormap(gray); axis off
        
        
        if OsampW(i)<.5
            % scale colormap based on range of nearest (larger) optimal JE
            subplot(2,nJ,i);
            clim([min(WsOld(:)),max(WsOld(:))]);
            
            [~,D1] = eig(Ws(1:NactW(i  ),1:NactW(i  )));
            [~,D2] = eig(Ws(1:NactW(i-1),1:NactW(i-1)));
            subplot(2,nJ,i+nJ);hold on;
            bar([1,2],[max(diag(D1)),max(diag(D2))]);
            plot([0,3],[1,1],'--k')
            xlim([0,3])
            ylim([.25,1.75])
            xticks([]);
            yticks([.5,1,1.5])
            yticklabels({'-','0','+'})
            set(gca,'fontsize',fs)
        else
            [~,D] = eig(Ws(1:NactW(i),1:NactW(i)));
            subplot(2,nJ,i+nJ);hold on;
            bar(1.5,max(diag(D)));
            plot([0,3],[1,1],'--k')
            xlim([0,3])
            ylim([.25,1.75])
            xticks([]);
            yticks([.5,1,1.5])
            yticklabels({'-','0','+'})
            set(gca,'fontsize',fs)
        end
        WsOld = Ws-min(Ws(:))-dW;
    end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(7,plotInds)
    %plot drift trajectories with stable/unstable regimes superimposed
    
    figure;set(gcf,'Position',[200 200 1600 800],'color','w');
    
    Jtmp = [JEsampW(find(JEsampW>JEsamp0,1,'last')),JEsamp0];
    Otmp = [1,0];
    for j=1:numel(Jtmp)
        subplot(2,3,j);hold on;
        [actNeurons,indsAct,~,~] = plotMultDriftTraj_regimes(Jtmp(j),JEmin,JEmax,N,tmax1,tmax2,psi0,Otmp(j),1,simVars,plotVars);
        
        nact = unique(sum(actNeurons));
        subplot(2,3,j+3);hold on;
        xlim([psi0(indsAct(1)),psi0(indsAct(end))])
        axis off
        for i=1:numel(indsAct)
            n = sum(actNeurons(:,i));
            ii = find(actNeurons(:,i));
            plot(psi0(indsAct(i)),th,'ok')
            if Otmp(j)
                plot(psi0(indsAct(i)),th(ii),'o','markerfacecolor','k','markeredgecolor','none')
            else
                if n>mean(nact)
                    plot(psi0(indsAct(i)),th(ii),'o','markerfacecolor',cO,'markeredgecolor','none')
                else
                    plot(psi0(indsAct(i)),th(ii),'o','markerfacecolor',cT,'markeredgecolor','none')
                end
            end
        end
        
    end
    subplot(2,3,3);hold on;
    plotSingleDriftTraj(Jtmp(2),JEmin,JEmax,N,tmax1,tmax2,simVars,plotVars);
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(8,plotInds)
    %plot drift rates and widths of stable and unstable regimes  
    figure;set(gcf,'Position',[200 200 700 400],'color','w');
    plotHeatmap_color(JEmin,JEmax,N,simVars,plotVars)
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(9,plotInds)
    %plot simplified energy in absence of velocity input
    figure;set(gcf,'Position',[200 200 800 600],'color','w');
    hold on;

    subplot(1,2,1);hold on;
    plotSimplifiedEnergy(JEsamp0,JEmin,JEmax,N,0,simVars,plotVars,0,1)
    xlabel('bump orientation')
    ylabel('energy')

    subplot(1,2,2);hold on;
    for i=1:numel(JEsampE)
        plotSimplifiedEnergy(JEsampE(i),JEmin,JEmax,N,0,simVars,plotVars,1./JEsampE(i),.1)
        ylbl{i} = ['1/',num2str(round(JEsampE(i),2))];
        plot([-pi/N,pi/N],[1./JEsampE(i),1./JEsampE(i)],'--k')
    end
    plot([-pi/N,0],[1./JEmax,1./JEmin],'--k')
    plot([0,pi/N ],[1./JEmin,1./JEmax],'--k')
    yticks(1./JEsampE)
    yticklabels(ylbl)
    xlabel('bump orientation')
    ylabel('energy')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(10,plotInds)
    % plot bump trajectory in different regimes, without input
    
    figure;set(gcf,'Position',[100 100 1400 200],'color','w');
    
    JEtmp = fliplr(JEsampN);
    for j=1:numel(JEtmp)
        subplot(1,5,j);hold on;
        plotMultDriftTraj_regimes(JEtmp(j),JEmin,JEmax,N,tmax1,tmax2,psi0,0,1,simVars,plotVars);
        
    end
        
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(11,plotInds) 
    %compare drift speed for analytics vs simulations
    
    vin = 0;
    t = 0:dt:30;
    
    figure;set(gcf,'Position',[200 200 700 400],'color','w');hold on
    
    JEint = 1./linspace(1./JEmax,1./JEmin,100);
    for i=1:numel(JEint)
        [~,~,tauU,tauS] = computeTimeConstants(JEint(i),JEmin,JEmax,tau);
        tDrift(i)    = computeDriftTime(tauU,tauS);
        driftRate(i) = computeNetDriftRate(tDrift(i),N);
    end
    plot(1./JEint,driftRate,'-k','linewidth',lw)
    xlim([1./JEmax,1./JEmin])
    ylim([0,.5])
    
    for i=1:numel(JEsim)
        [~,~,tauU,tauS] = computeTimeConstants(JEsim(i),JEmin,JEmax,tau);
        [dpsiU,dpsiS] = computeDelNact(tauU,tauS,N);
        psi0 = pi/N-dpsiU/exp(1);
        psiF = dpsiS/exp(1);
        
        [psi,~] = simPsi(JEsim(i),vin,th,A,C0,psiSamp,tau,N,t,tV,del,psi0);
        psi = mod(psi,2*pi);
        ii = find(psi<psiF,1,'first');
        tDriftSim(i)    = t(ii);
        driftRateSim(i) = computeNetDriftRate(tDriftSim(i),N);
    end
    plot(1./JEsim,driftRateSim,'ok','markersize',ms)
    xlabel('excitation 1/J_E')
    ylabel('net drift speed')
    [ticks,labels] = getJElabels(JEmin,JEmax,7);
    xticks(ticks)
    xticklabels(labels);
    set(gca,'fontsize',fs)
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(12,plotInds) 
    %plot simplified energy in presence of velocity input

    figure;set(gcf,'Position',[200 200 800 600],'color','w');
    hold on;

    subplot(1,2,1);hold on
    vin = [.15,.3,.45];
    for i=1:numel(vin)
        plotSimplifiedEnergy(JEsamp0,JEmin,JEmax,N,vin(i),simVars,plotVars,0,1)
    end
    xlabel('bump orientation')
    ylabel('energy')

    subplot(1,2,2);hold on;
    vin = .125;
    for i=1:numel(JEsampE)
        plotSimplifiedEnergy(JEsampE(i),JEmin,JEmax,N,vin,simVars,plotVars,1./JEsampE(i),.1)
        plot([-pi/N,pi/N],[1./JEsampE(i),1./JEsampE(i)],'--k')
        ylbl{i} = ['1/',num2str(round(JEsampE(i),2))];
    end
    plot([-pi/N,0],[1./JEmax,1./JEmin],'--k')
    plot([0,pi/N ],[1./JEmin,1./JEmax],'--k')

    JEint = 1./linspace(1./JEmax,1./JEmin,100);
    [lambdaU,lambdaS,~,~] = computeTimeConstants(JEint,JEmin,JEmax,tau);
    psiU = computeUnstableFPs(vin,lambdaU,N);
    psiS = computeStableFPs(  vin,lambdaS);
    plot(psiU,1./JEint,'--k','linewidth',lw)
    plot(psiS,1./JEint,'-k','linewidth',lw)

    yticks(1./JEsampE)
    yticklabels(ylbl)
    xlabel('bump orientation')
    ylabel('energy')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(13,plotInds)
    %plot bump trajectory in different regimes, with velocity input 

    figure;set(gcf,'Position',[200 200 1400 200],'color','w');
    
    vin = mdesired_bumpTraj_varJE;
    JEint = fliplr(JEsampN);
    psiInit = 0;
    for i=1:numel(JEint)
        
        subplot(1,5,i);hold on;
        plotMultDrivenTrajAboveThresh_regimes(vin,JEint(i),JEmin,JEmax,N,tmax2,tmax2,0,1,simVars,plotVars,psiInit)
         
    end
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(14,plotInds)
    %plot threshold velocity
    
    figure;set(gcf,'Position',[200 200 700 400],'color','w');hold on;

    %plot threshold velocity (analytic & simulation)
    JEint = 1./linspace(1./JEmax,1./JEmin,1000);
    [~,~,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
    [~,dpsiS] = computeDelNact(tauU,tauS,N);
    vmin = computeVmin(tauS,dpsiS);
    plot(1./JEint,vmin,'-k','linewidth',lw)
    
    [ticks,labels] = getJElabels(JEmin,JEmax,7);
    xticks(ticks);
    xticklabels(labels);
    xlim([1./JEmax,1./JEmin])  
    ylabel('input velocity')
    set(gca,'fontsize',fs,'tickdir','out')
    
    for i=1:numel(JEsim)
        vminSim(i) = simVmin(JEsim(i),JEmin,JEmax,N,tmax2,simVars);
    end
    plot(1./JEsim,vminSim,'ok','markersize',ms)
    ylim([0,.8])
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(15,plotInds)
    %plot linearity of itnegration
    figure;set(gcf,'Position',[200 200 1400 400],'color','w');
     
    %plot FPs above threshold
    subplot(1,3,1);hold on;
    plotHeatmap(JEmin,JEmax,N,simVars,plotVars)
    JEint = 1./linspace(1./JEmax,1./JEmin,1000);
    [~,~,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
    [~,dpsiS] = computeDelNact(tauU,tauS,N);
    vmin = computeVmin(tauS,dpsiS);
    
    %plot linearity of integration (analytic & simulation)
    nV = numel(mdesired_linearity);
    cmap = gray(nV+2);
    cmap = flipud(cmap(1:end-2,:));
    
    subplot(1,2,1);hold on;
    for i=1:numel(JEsim)
        linS(i,1:nV) = simLinInt(mdesired_linearity,JEsim(i),JEmin,JEmax,N,tmax1,tmax2,simVars);
    end
    
    for i=1:nV
        linA = computeLinearity(mdesired_linearity(i),vmin);
        plot(1./JEint,linA,'-','linewidth',lw,'color',cmap(i,:));
        
        plot(1./JEsim,linS(:,i),'o','markersize',ms,'markeredgecolor',cmap(i,:),'markerfacecolor','none')
    end
    [ticks,labels] = getJElabels(JEmin,JEmax,7);
    xticks(ticks);
    xticklabels(labels);
    xlim([1./JEmax,1./JEmin])  
    ylabel('linearity')
    set(gca,'fontsize',fs,'tickdir','out')
    
    %plot velocity colormap
    subplot(1,2,2);hold on;
    plotLineColors(mdesired_linearity,plotVars)
   
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(16,plotInds)
    %plot robustness to parameter variations
    
    nSamples = 6:2:20;
    var = 'vmin';
    
    if strcmp(var,'drift')
        xlm = [0,1.2];
        fun = @log;
        cax = [-5,1];
        thresh = 1e-3;
    elseif strcmp(var,'vmin')
        xlm = [0,1.5];
        fun = @log;
        cax = [-7,1];
        thresh = 1e-3;
    elseif strcmp(var,'lin')
        xlm = [0,1];
        fun = @(x) x;
        cax = [0,1];
        thresh = 1e-4;
        thresh = 1-thresh;
    else
        error('unrecognized input variable')
    end
    
    fig1 = figure();
    set(fig1,'Position',[200 200 1800 800],'color','w');
    
    plotNdepHeatmap(var,nSamples,simVars,plotVars,fig1,[2,3,1]);
    yticklabels({})
    ylabel('')

    plotNdepThresh(var,thresh,nSamples,simVars,plotVars,fig1,[2,3,2]);
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(17,plotInds)
    %plot noise robustness
    
    noiseLevels = A.*[1/6, 2/6, 3/6, 4/6, 5/6];

    nN = 6;
    %uncomment to generate and save simulations for a network of nN neurons
    %nN = 6:2:20;
    %for i=1:numel(nN)
    %   genNoiseVar_varJE(nN(i),noiseLevels(1),simVars);
    %end
    plotNoiseVar_varJE(nN,noiseLevels(1),plotVars);

   
    %uncomment to generate and save simulations for different noise levels
    %for i=1:numel(noiseLevels)
    %   genNoiseVar_varN(JEopt(2),noiseLevels(i),simVars);
    %end
    plotNoiseVar_varN(JEopt(2),noiseLevels,plotVars);
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(18,plotInds)
    %plot phase diagram

    psiFP = [th(1), 1/2*(th(1)+th(2))]; % fixed point values for psi
    JEsampD = 1./linspace(1/20,1/2.03,100); % JE sampling (D for dense)
    JIbif = zeros(2,length(JEsampD));

    for je = 1:length(JEsampD) % computing JIbif

        % compute fixed point values for thc
        thcFP = zeros(size(psiFP));
        thcinit = pi/2;
        [thcFP(1),fval,flag,~] = fzero(@(thc)([0 JEsampD(je) 0]*calculateFunctions([psiFP(1) thc],N)' - 1),thcinit); % compute value of thcstar
        if flag < 0
            disp(['JE = ' num2str(JEsampD(je)) ', 1: Something went wrong, flag for root finder: ' num2str(flag)])
        end
        [thcFP(2),fval,flag,~] = fzero(@(thc)([0 JEsampD(je) 0]*calculateFunctions([psiFP(2) thc],N)' - 1),thcinit);
        if flag < 0
            disp(['JE = ' num2str(JEsampD(je)) ', 2: Something went wrong, flag for root finder: ' num2str(flag)])
        end

        [f,~] = calculateFunctions([psiFP', thcFP'],N);
        f0 = f(:,1);

        JIbif(:,je) = -cos(thcFP')./f0;

    end

    figure
    set(gcf,'Position',[200 200 500 500],'color','w')
    plot(min(JIbif),1./JEsampD,'k','linewidth',lw)
    ylim([0, 0.6])
    xlim([-11,1.5])

    line([1,1],[0.6,1/2], 'color','k','linewidth',lw)
    line([min(min(JIbif)),1],[1/2,1/2], 'color','k','linewidth',lw)

    for je = 1:length(JEopt)
        thcSamp = zeros(size(psiSamp));
        for p = 1:length(psiSamp)
            thcinit = pi/2;
            [thcSamp(p),fval,flag,~] = fzero(@(thc)([0 JEopt(je) 0]*calculateFunctions([psiSamp(p) thc],N)' - 1),thcinit);
            if flag < 0
                disp(['JE = ' num2str(JEsampD(je)) ', p = ' num2str(psiSamp(p)) ', something went wrong, flag for root finder: ' num2str(flag)])
            end

            [f,~] = calculateFunctions([psiSamp', thcSamp'],N);
            f0 = f(:,1);

            JI_bound_opt = min(-cos(thcSamp')./f0);

        end
        line([min(min(JIbif)),JI_bound_opt],[1/JEopt(je),1/JEopt(je)],'color','k','linestyle','--','linewidth',lw)
    end

    ylabel('local excitation 1/J_E')
    xlabel('broad inhibition J_I')
    title('stability of population profile')
    set(gca,'fontsize',fs,'tickdir','out')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(19,plotInds)
    %plot fixed point conditions

    JEinv = linspace(1/JEopt(end),1/JEopt(1),10);

    nDense = 1000; % number of samples for dense sampling
    thcSampD = linspace(pi/N, (N-1)*pi/N, nDense+1); % half-width samples
    thcSampD(end) = []; thcSampD = thcSampD';
    psiSampD = linspace(0,2*pi,nDense+1); % psi samples (D for dense sampling)
    psiSampD(end) = []; psiSampD = psiSampD';

    % space allocation
    feven = zeros(length(thcSampD),length(psiSampD));
    fodd = zeros(length(thcSampD),length(psiSampD));

    for p = 1:length(psiSampD)
        [f,~] = calculateFunctions([psiSampD(p)*ones(size(thcSampD)),thcSampD],N);
        feven(:,p) = f(:,2);
        fodd(:,p) = f(:,3);
    end

    cmap = magma(length(JEinv)+1); % colormap for feven plot
    cmap = cmap(1:end-1,:);

    figure
    set(gcf,'Position',[200 200 1000 500],'color','w')
    clear ax

    ax(1) = subplot(1,2,1);
    imagesc(psiSampD,2*thcSampD,fodd);
    colormap(ax(1),flipud(redblueu))
    cb = colorbar('northoutside');
    set(gca,'Ydir','normal')
    xlabel('bump orientation \psi (rad)')
    ylabel('bump width w')
    cb.Label.String = 'f_{odd}';
    set(gca,'fontsize',fs,'tickdir','out')

    ax(2) = subplot(1,2,2);
    contour(psiSampD, 2*thcSampD, feven, JEinv,'linewidth',lw);
    cb = colorbar('northoutside');
    cb.Ticks = JEinv(1:3:end);
    cb.Label.String = '1/J_E';
    colormap(ax(2),cmap)
    yticklabels('')
    xlabel('bump orientation \psi (rad)')
    set(gca,'fontsize',fs,'tickdir','out')

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(20,plotInds)
    %plot linearization

    psiFP = [th(1), 1/2*(th(1)+th(2))]; % fixed point values for psi
    JEsampD = 1./linspace(1/20,1/2.03,500); % JE sampling (D for dense)
    JIsamp = linspace(-11,1);

    % space allocation
    lamPsi = zeros(length(JEsampD),length(psiFP));
    JIBound = zeros(length(JEsampD),length(psiFP));
    lamPlus = zeros(length(JEsampD),length(JIsamp),length(psiFP));
    lamMinus = zeros(length(JEsampD),length(JIsamp),length(psiFP));

    for p = 1:length(psiFP)
        psi = psiFP(p);
        for je = 1:length(JEsampD)
            JE = JEsampD(je);

            % initial value for thetac
            thcinit = pi/2;
            [thc,~,flag,~] = fzero(@(thc)([0 JE 0]*calculateFunctions([psi thc],N)' - 1),thcinit);
            if flag < 0
                disp(['Something went wrong, flag for root finder: ' num2str(flag)])
            end

            [f,actInd] = calculateFunctions([psi thc], N);
            Nact = sum(actInd);
            JIBound(je,p) = -cos(thc)/f(1);

            actInd = logical(actInd);

            lamPsi(je,p) = -1+JE/N*sum(sin(th(actInd) - psi).^2);
            cosSum = sum(cos(th(actInd) - psi));

            for ji = 1:length(JIsamp)
                JI = JIsamp(ji);
                a = -C0/(cos(thc) + JI*f(1));

                discrim = ((JI - JE)*Nact/N + lamPsi(je,p) + 1)^2 + 4*JI*JE*cosSum^2/N^2;
                linTerm = (JI + JE)*Nact/N - 3 - lamPsi(je,p);

                lamPlus(je,ji,p) = 1/2*(linTerm + sqrt(discrim));
                lamMinus(je,ji,p) = 1/2*(linTerm - sqrt(discrim));

            end
        end
    end

    figure
    set(gcf,'Position',[200 200 1200 800],'color','w')
    clear ax

    ax(1) = subplot(2,3,1);
    plot(lamPsi(:,1),1./JEsampD,'k','linewidth',lw)
    ylabel({'\psi^* = \theta_j','','local excitation, 1/J_E'})
    xticklabels('')
    
    ax(2) = subplot(2,3,4);
    plot(lamPsi(:,2),1./JEsampD,'k','linewidth',lw)
    ylabel({'\psi^* = 1/2(\theta_j + \theta_{j+1})','','local excitation, 1/J_E'})
    xlabel('\lambda_\psi')

    for a = 1:length(ax)
        line(ax(a),[0,0],1./[JEsampD(1),JEsampD(end)],'color','k','linestyle','--')
        for i = 1:length(JEopt)
            line(ax(a),[-1.5,2.5],1./JEopt(i)*[1,1],'color','k','linestyle',':')
        end
        xlim(ax(a),[-1.5,2.5])
        ylim(ax(a), 1./[JEsampD(1),JEsampD(end)])
        set(ax(a),'fontsize',fs,'tickdir','out')
    end
    
    clear ax

    ax(1) = subplot(2,3,2); hold on;
    imagesc(JIsamp,1./JEsampD,lamPlus(:,:,1))
    colormap(ax(1),flipud(redblueu([-4.5,1])))
    clim([-1,4.5])
    cb = colorbar('northoutside');
    cb.Label.String = '\lambda_+';
    yticklabels('')
    xticklabels('')

    ax(2) = subplot(2,3,5); hold on;
    imagesc(JIsamp,1./JEsampD,lamPlus(:,:,2))
    colormap(ax(2),flipud(redblueu([-4.5,1])))
    clim([-1,4.5])
    xlabel('broad inhibition J_I')
    yticklabels('')

    ax(3) = subplot(2,3,3); hold on;
    imagesc(JIsamp,1./JEsampD,lamMinus(:,:,1))
    colormap(ax(3),flipud(redblueu([0,10.5])))
    clim([-10.5,0])
    cb = colorbar('northoutside');
    cb.Label.String = '\lambda_-';
    yticklabels('')
    xticklabels('')

    ax(4) = subplot(2,3,6); hold on; 
    imagesc(JIsamp,1./JEsampD,lamMinus(:,:,2))
    colormap(ax(4),flipud(redblueu([0,10.5])))
    clim([-10.5,0])
    xlabel('broad inhibition J_I')
    yticklabels('')

    for a = 1:length(ax)
        set(ax(a),'ydir','normal')
        plot(ax(a),JIBound(:,1),1./JEsampD,'k','linewidth',lw)
        for i = 1:length(JEopt)
            line(ax(a),[JIsamp(1),JIsamp(end)],1./JEopt(i)*[1,1],'color','k','linestyle',':')
        end
        ylim(ax(a),1./[max(JEsampD),min(JEsampD)])
        xlim(ax(a),[min(JIsamp),max(JIsamp)])
        set(ax(a),'fontsize',fs,'tickdir','out')
    end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(21,plotInds)
    %plot eigenvalues of Hessian
    
    th0 = th';
    JEopt = computeOptJE(N);
    psiSamp0 = linspace(th(1),th(2),400);
    nJ = numel(JEopt);
    
    
    figure
    set(gcf,'Position',[200 200 1200 800],'color','w')
    hold on
        
    for i=1:4*nJ
        subplot(4,nJ,i);hold on
        for k=[0,pi/N,2*pi/N]
            plot([k,k],[-1,1],'--r')
        end
    end
        
        
    for j=1:numel(JEopt)
        JE = JEopt(j);
        [~,~,JI] = getConnectivity(N,JE,A,C0,psiSamp0,del);


        for p = 1:length(psiSamp0)

            psi = psiSamp0(p);

            thcinit = pi/2;
            [thc,~,flag,~] = fzero(@(thc)([0 JE 0]*calculateFunctions0([psi thc],th0,N) - 1),thcinit);
            if flag < 0
                disp(['Something went wrong, flag for root finder: ' num2str(flag)])
            end

            [f,~] = calculateFunctions0([psi thc],th0,N);
            f0   = f(1);
            rho = -C0/(2*(cos(thc) + JI*f0));

            H = calculateSecondPartials(rho,thc,psi,JI,JE,C0,th0,N);

            [V,D] = eig(H);
            lam = diag(D);
            
            rhoDir = find(abs(V(1,:)) == max(abs(V(1,:))));
            thcDir = find(abs(V(2,:)) == max(abs(V(2,:))));
            psiDir = find(abs(V(3,:)) == max(abs(V(3,:))));


            if numel(psiDir)>1
                psiDir(psiDir==rhoDir) = [];
                psiDir(psiDir==thcDir) = [];
            end

            lam1(j,p) = lam(rhoDir);
            lam2(j,p) = lam(thcDir);
            lam3(j,p) = lam(psiDir);

            v1rho(j,p) = V(1,rhoDir);
            v1thc(j,p) = V(2,rhoDir);
            v1psi(j,p) = V(3,rhoDir);

            v2rho(j,p) = V(1,thcDir);
            v2thc(j,p) = V(2,thcDir);
            v2psi(j,p) = V(3,thcDir);

            v3rho(j,p) = V(1,psiDir);
            v3thc(j,p) = V(2,psiDir);
            v3psi(j,p) = V(3,psiDir);
        end
        
        
        if j<nJ
            inds = 2:numel(psiSamp0)-1;
            subplot(4,nJ,j);hold on;plot(psiSamp0(inds),lam3(j,inds),'-k','linewidth',lw)
            ylim([-1,1]); 
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,nJ+j); hold on; plot(psiSamp0(inds),abs(v3psi(j,inds)),'-k','linewidth',lw)
            ylim([.94,1]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,2*nJ+j);hold on;plot(psiSamp0(inds),abs(v3thc(j,inds)),'-k','linewidth',lw)
            ylim([0,.4]);yticks([0,.2,.4]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,3*nJ+j);hold on;plot(psiSamp0(inds),abs(v3rho(j,inds)),'-k','linewidth',lw)
            ylim([0,.3]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({'0','\pi/N','2\pi/N'})
            set(gca,'fontsize',fs,'tickdir','out')
            


        else
            for i=2:numel(psiSamp0)-1
                A=null([v1rho(j,i);v1thc(j,i);v1psi(j,i)]');
                vv2rho(j,i) = A(1,1);
                vv2thc(j,i) = A(2,1);
                vv2psi(j,i) = A(3,1);
                vv3rho(j,i) = A(1,2);
                vv3thc(j,i) = A(2,2);
                vv3psi(j,i) = A(3,2);
            end

            
            inds = 2:numel(psiSamp0)-1;
            subplot(4,nJ,j);hold on;plot(psiSamp0(inds),lam3(j,inds),'-k','linewidth',lw)
            ylim([-1,1]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,nJ+j);hold on;  plot(psiSamp0(inds),abs(vv3psi(j,inds)),'-k','linewidth',lw)
            ylim([.94,1]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,2*nJ+j);hold on;plot(psiSamp0(inds),abs(vv3thc(j,inds)),'-k','linewidth',lw)
            ylim([0,.4]);yticks([0,.2,.4]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            subplot(4,nJ,3*nJ+j);hold on;plot(psiSamp0(inds),abs(vv3rho(j,inds)),'-k','linewidth',lw)
            ylim([0,.3]);
            xlim([0,2*pi/N]);xticks([0,pi/N,2*pi/N]);xticklabels({'0','\pi/N','2\pi/N'})
            set(gca,'fontsize',fs,'tickdir','out') 

        end
        
    end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(22,plotInds)
    %plot leading eigenvalues, no input
    figure;set(gcf,'Position',[200 200 1200 800],'color','w');
    
    nSamples = [6,8,10];
    nS = numel(nSamples);
    vin = 0;
    for i=1:nS
        [JEopt,~,Nact] = computeOptJE(nSamples(i));
        for j=1:numel(JEopt)-1
            JE = 1./linspace(1./JEopt(j+1),1./JEopt(j),101);
            [evU,evS] = computeEvals(JE,JEopt(j),JEopt(j+1),tau);
            subplot(1,nS,i);hold on;
            plot(evU,1./JE,'-','linewidth',lw,'color',cO)
            plot(evS,1./JE,'-','linewidth',lw,'color',cT)
            plot([-10,30],[1./JEopt(j),1./JEopt(j)],'--','color',cR)
            
            jsamps = 1:5:numel(JE);
            for k=1:numel(jsamps)
                [Ws,Wa,~] = getConnectivity(nSamples(i),JE(jsamps(k)),A,C0,psiSamp,del,1);
                W = getFullConnectivity(Ws,Wa,vin,tau);              
                [evUn(k),evSn(k)] = computeEvalsNum(W,Nact(j));
            end
            plot(evUn,1./JE(jsamps),'o','markeredgecolor',cO,'markerfacecolor','none','markersize',ms/2)
            plot(evSn,1./JE(jsamps),'o','markeredgecolor',cT,'markerfacecolor','none','markersize',ms/2)
        end
        plot([-10,30],[1./JEopt(j+1),1./JEopt(j+1)],'--','color',cR)
        subplot(1,nS,i);
        plot([0,0],[0,.5],'--k')
        title(['N = ',num2str(nSamples(i))])
        set(gca,'fontsize',fs,'tickdir','out')
        xlabel('leading eigenvalues')
        if i>1
            yticklabels({})
        else
            ylabel('excitation 1/J_E')
        end
    end
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(23,plotInds)
    %plot leading eigenvalues, small input
    
    figure;set(gcf,'Position',[200 200 400 600],'color','w');hold on;
    
    vin = 0:.1:1;
    nV = numel(vin);
    cmap = parula(nV+2);
    [JEopt,~,Nact] = computeOptJE(N);
    
    cmapU1 = buildColormap(cO,cW,(nV+1)/2+2);cmapU1(1:2,:) = [];   
    cmapU2 = buildColormap(cK,cO,(nV-1)/2);
    cmapU = [cmapU1;cmapU2];
    
    cmapS1 = buildColormap(cT,cW,(nV+1)/2+2);cmapS1(1:2,:) = [];   
    cmapS2 = buildColormap(cK,cT,(nV-1)/2);
    cmapS = [cmapS1;cmapS2];
    
    for i=1:nV
        for j=1:numel(JEopt)-1
            
            JE = 1./linspace(1./JEopt(j+1),1./JEopt(j),51);
            
            for k=1:numel(JE)
                [Ws,Wa,~] = getConnectivity(N,JE(k),A,C0,psiSamp,del,0);
                W = getFullConnectivity(Ws,Wa,vin(i),tau);
                
                [evUn(k),evSn(k)] = computeEvalsNum(W,Nact(j));
            end
            plot(evUn,1./JE,'-','linewidth',lw,'color',cmapU(i,:))
            plot(evSn,1./JE,'-','linewidth',lw,'color',cmapS(i,:))
        end
    end
    for j=1:numel(JEopt)	
        plot([-10,30],[1./JEopt(j),1./JEopt(j)],'--','color',cR)
    end
    plot([0,0],[0,.5],'--k')
    xlim([-10,20])
    ylim([1./max(JEopt),1./min(JEopt)])
    set(gca,'fontsize',fs,'tickdir','out')
    ylabel('excitation 1/J_E')
    xlabel('leading eigenvalues')
    
    cmap = buildColormap(cK,cW,100);
    cmap = [flipud(cmap);cmap];
    
    figure;set(gcf,'Position',[200 200 800 800],'color','w');hold on;
    vin = linspace(0,1,51);
    nV = numel(vin);
    nJ = numel(JEopt)-1;
    
    for j=1:nJ
        JE = 1./linspace(1./JEopt(j+1),1./JEopt(j),51);
        Js = 1./linspace(1./JEopt(j+1),1./JEopt(j),7);
        Js = Js([2,4,6]);
        for k=1:numel(JE)
            for i=1:nV            
                [Ws,Wa,~] = getConnectivity(N,JE(k),A,C0,psiSamp,del,0);
                W = getFullConnectivity(Ws,Wa,vin(i),tau);
                [evUn(k,i),evSn(k,i)] = computeEvalsNum(W,Nact(j));
            end
            ps(k,1:4) = polyfit(vin(2:end),(evSn(k,2:end)-evSn(k,1)),3);
            pu(k,1:4) = polyfit(vin(2:end),(evUn(k,2:end)-evUn(k,1)),3);
        end
        pall = [ps;pu];
        cmax = max(abs(pall(:)));

        subplot(6,4,12*(j-1)+1);imagesc(1:4,1./JE,fliplr(ps));colormap(cmap);clim([-cmax,cmax])
        set(gca,'fontsize',fs,'tickdir','out','ydir','normal')
 
        subplot(6,4,12*(j-1)+2);imagesc(1:4,1./JE,fliplr(pu));colormap(cmap);clim([-cmax,cmax])
        set(gca,'fontsize',fs,'tickdir','out','ydir','normal')
        
        vSmin = [];
        vUmin = [];
        for k = 1:numel(Js)
            [~,kk] = min(abs(Js(k)-JE));
            vS = evSn(kk,2:end)-evSn(kk,1);
            vU = evUn(kk,2:end)-evUn(kk,1);
            
            vSfit = ps(kk,1).*vin(2:end).^3+ps(kk,2).*vin(2:end).^2;
            vUfit = pu(kk,1).*vin(2:end).^3+pu(kk,2).*vin(2:end).^2;
            
            vSmin = [vSmin,min(vS)];
            vUmin = [vUmin,min(vU)];
            
            subplot(6,4,12*(j-1)+3+4*(k-1));hold on;
            plot(vin(2:end),vS,     'linewidth',lw,'color',cT)
            plot(vin(2:end),vSfit,'--','linewidth',lw,'color',.7*cT)
            
            subplot(6,4,12*(j-1)+4+4*(k-1));hold on;
            plot(vin(2:end),vU,        'linewidth',lw,'color',cO)
            plot(vin(2:end),vUfit,'--','linewidth',lw,'color',.7*cO)
        end
        
        if max(abs(vSmin))>.1
            fS = 20;
        else
            fS = 200;
        end
        if max(abs(vUmin))>1
            fU = 2;
        else
            fU = 20;
        end
        for k = 1:numel(Js)
            subplot(6,4,12*(j-1)+3+4*(k-1))
            ylim([floor(fS*min(vSmin))./fS,0])
            xticklabels({})
            set(gca,'fontsize',fs,'tickdir','out')
            
            
            subplot(6,4,12*(j-1)+4+4*(k-1))
            ylim([floor(fU*min(vUmin))./fU,0])
            set(gca,'fontsize',fs,'tickdir','out')
            xticklabels({})
        end
        
        
    end
    
   
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(24,plotInds)
    %plot drift rates and widths of stable and unstable regimes  
    figure;set(gcf,'Position',[200 200 700 400],'color','w');
    plotHeatmap(JEmin,JEmax,N,simVars,plotVars)
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(25,plotInds)
    %plot fixed points in different regimes for varying input velocity
    
    figure;set(gcf,'Position',[200 200 2000 600],'color','w');
    nV = numel(mdesired_bumpTraj);
    for i=1:nV
        subplot(2,nV,i);
        plotHeatmap(JEmin,JEmax,N,simVars,plotVars)
        
        hold on;
        %plot FPs for all JE
        JEint = 1./linspace(1./JEmax,1./JEmin,100);
        [lambdaU,lambdaS,~,~] = computeTimeConstants(JEint,JEmin,JEmax,tau);
        psiU = computeUnstableFPs(mdesired_bumpTraj(i),lambdaU,N);
        psiS = computeStableFPs(  mdesired_bumpTraj(i),lambdaS);
        plot(psiU,1./JEint,'--k','linewidth',lw)
        plot(psiS,1./JEint,'-k','linewidth',lw)
        xlabel('')
        xticklabels({})
        if i>1
            yticklabels({})
            ylabel('')
        end
        
        %plot FPs for single JE
        [lambdaU,lambdaS,~,~] = computeTimeConstants(JEsamp0,JEmin,JEmax,tau);
        psiU = computeUnstableFPs(mdesired_bumpTraj(i),lambdaU,N);
        psiS = computeStableFPs(  mdesired_bumpTraj(i),lambdaS);
        plot([0,pi/N],[1./JEsamp0,1./JEsamp0],'--k')
        plot(psiS,1./JEsamp0,'s','markersize',ms,'markerfacecolor','k','markeredgecolor','none')
        plot(psiU,1./JEsamp0,'s','markersize',ms,'markerfacecolor','w','markeredgecolor','k')

        subplot(2,nV,i+nV);hold on;
        plotDoubleDrivenTraj(mdesired_bumpTraj(i),JEsamp0,JEmin,JEmax,N,tmax1,tmax2,simVars,plotVars)
        if i>1
            yticklabels({})
            ylabel('')
        end
    end
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(26,plotInds)
    %plot fixed points in different regimes for varying input velocity
    %below threshold
    
    figure;set(gcf,'Position',[200 200 700 400],'color','w');
    
    %plot FPs below threshold
    alphas = linspace(0,1,numel(mdesired_FPbelowThresh)+2);
    alphas = alphas(3:end);
    hold on;
    %plotHeatmap(JEmin,JEmax,N,simVars,plotVars)
    JEint = 1./linspace(1./JEmax,1./JEmin,1000);
    [lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
    [~,dpsiS] = computeDelNact(tauU,tauS,N);
    for i=1:numel(mdesired_FPbelowThresh)
        psiU = computeUnstableFPs(mdesired_FPbelowThresh(i),lambdaU,N);
        psiS = computeStableFPs(  mdesired_FPbelowThresh(i),lambdaS);
        ii = find(psiS<dpsiS);
        jj = find(psiU>dpsiS);
        plot(1./JEint(ii),psiS(ii),'-','color',[cT,alphas(i)],'linewidth',lw)
        plot(1./JEint(jj),psiU(jj),'--','color',[cO,alphas(i)],'linewidth',lw)  
    end
    ylim([0,pi/N])
    yticks([0,pi/(2*N),pi/N])
    yticklabels({'0','\pi/2N','\pi/N'})
    xlim([1./JEmax,1./JEmin])
    [ticks,labels] = getJElabels(JEmin,JEmax,7);
    xticks(ticks);
    xticklabels(labels);
    set(gca, 'fontsize',fs)
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
if ismember(27,plotInds)
    %plot fixed points in different regimes for varying input velocity
    %below threshold

    figure;set(gcf,'Position',[200 200 700 400],'color','w');
     
    %plot FPs above threshold
    alphas = linspace(0,1,numel(mdesired_FPaboveThresh)+2);
    alphas = alphas(3:end);
    hold on;
    %plotHeatmap(JEmin,JEmax,N,simVars,plotVars)
    JEint = 1./linspace(1./JEmax,1./JEmin,1000);
    [lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
    [~,dpsiS] = computeDelNact(tauU,tauS,N);
    vmin = computeVmin(tauS,dpsiS);
    for i=1:numel(mdesired_FPaboveThresh)
        
        psiU = computeUnstableFPs(mdesired_FPaboveThresh(i),lambdaU,N);
        psiS = computeStableFPs(  mdesired_FPaboveThresh(i),lambdaS);
        ii = find(psiS>dpsiS);
        jj = find(psiU<dpsiS);
        plot(1./JEint(ii),psiS(ii),'-','color',[cT,alphas(i)],'linewidth',lw)
        plot(1./JEint(jj),psiU(jj),'--','color',[cO,alphas(i)],'linewidth',lw)  
    end
    
    ylim([0,pi/N])
    yticks([0,pi/(2*N),pi/N])
    yticklabels({'0','\pi/2N','\pi/N'})
    xlim([1./JEmax,1./JEmin])
    [ticks,labels] = getJElabels(JEmin,JEmax,7);
    xticks(ticks);
    xticklabels(labels);
    set(gca, 'fontsize',fs)
   
end
%-------------------------------------------------------------------------%



end
%-------------------------------------------------------------------------%






%------------------- SETTING PLOT FEATURES -------------------------------%

function plotLineColors(vec,plotVars)

lw = plotVars.lw;
fs = plotVars.fs;

nV = numel(vec);
colors = gray(nV+2);
colors = flipud(colors(1:end-2,:));

% plot velocity colors:
dv = vec(2)-vec(1);

for i=1:numel(vec)
    plot([vec(i),vec(i)],[1,2],'linewidth',lw,'color',colors(i,:));
end
xlabel('input velocity')
xlim([vec(1)-dv,vec(end)+dv])
xticks(vec)
ylim([0,3])
yticks([])
set(gca,'fontsize',fs,'TickDir','out');

end

function [ticks,labels] = getJElabels(Jmin,Jmax,n)
ticks = linspace(1./Jmax,1./Jmin,n);
for i=1:n
    labels{i} = ['1/',num2str(round(1./ticks(i),2))];
end
end

function [cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,np)
[~,lSmax,~,~] = computeTimeConstants(JEmin,JEmin,JEmax,tau);
[lUmax,~,~,~] = computeTimeConstants(JEmax,JEmin,JEmax,tau);
[lU,lS,~,~  ] = computeTimeConstants(JE,   JEmin,JEmax,tau);

lSall = linspace(0,lSmax,np);
[~,indS] = min((lS-lSall).^2);
cmapT = buildColormap(cT,cW,np);
cS = cmapT(indS,:);

lUall = linspace(0,lUmax,np);
[~,indU] = min((lU-lUall).^2);
cmapO = buildColormap(cO,cW,np);
cU = cmapO(indU,:);

%NEW: single B&W colormap
lall = linspace(0,max(lUmax,abs(lSmax)),np);
[~,indU] = min((lU-lall).^2);
[~,indS] = min((abs(lS)-lall).^2);
cmap = buildColormap([0,0,0],cW,np);
cU = cmap(indU,:);
cS = cmap(indS,:);
end

function cmap = buildColormap(c1,c2,np)
% linearly interpolates between c1 and c2
% c1,c2: RGB values (c2 lighter (larger) than c1)
% np: number of points to sample colormap

cmap = [];
for i=1:np
    cmap = [cmap;c1+i*(c2-c1)./np];
end
cmap = flipud(cmap);

end


%--------------------- PLOTTING FUNCTIONS --------------------------------%

function plotSimplifiedEnergy(JE,JEmin,JEmax,N,vin,simVars,plotVars,offset,scale)

if nargin<9
    scale = 1;
end

tau = simVars.tau;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;

[lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
psiU = computeUnstableFPs(vin,lambdaU,N);
psiS = computeStableFPs(  vin,lambdaS);
[dpsiU,dpsiS] = computeDelNact(tauU,tauS,N);

psiU0 = computeUnstableFPs(0,lambdaU,N);
psiS0 = computeStableFPs(  0,lambdaS);

psiax = linspace(-psiU0-dpsiU,psiU0+dpsiU,1000);
iS = find(psiax>=psiS0-dpsiS & psiax<=psiS0+dpsiS);
iU1 = find(psiax<=psiS0-dpsiS);
iU2 = find(psiax>=psiS0+dpsiS);


psiU1 = (-2*psiU0+psiU);
psiU2 = psiU;
dpsi = 2*pi/N;

C2 = dpsi*(lambdaU*dpsiU/2 - vin + vin.^2./(lambdaU.*dpsiU.*2));
C1 = (dpsiU.^2*lambdaU-dpsiS.^2*lambdaS) + vin.*dpsi + vin.^2.*(1./lambdaU-1./lambdaS);

ES  = -lambdaS.*(psiax-psiS).^2;
EU1 = -lambdaU.*(psiax-psiU1).^2+C1;
EU2 = -lambdaU.*(psiax-psiU2).^2+C2;
vshift = ES(iS(end));
ES  = scale*(ES  - vshift);
EU1 = scale*(EU1 - vshift);
EU2 = scale*(EU2 - vshift);

%vertically shift so that that transition between stable and unstable
%regimes is at energy=0

plot(psiax(iS), ES( iS )+offset,'color',cT,'linewidth',lw);
plot(psiax(iU1),EU1(iU1)+offset,'color',cO,'linewidth',lw);
plot(psiax(iU2),EU2(iU2)+offset,'color',cO,'linewidth',lw);
set(gca,'fontsize',fs)


[~,inds ] = min(abs(psiS -psiax));
[~,indu1] = min(abs(psiU1-psiax));
[~,indu2] = min(abs(psiU2-psiax));
plot( psiS, ES(inds)  +offset,'s','markersize',12,'markerfacecolor',cT,'markeredgecolor','none')  
plot(psiU1, EU1(indu1)+offset,'s','markersize',12,'markerfacecolor',cO,'markeredgecolor','none')
plot(psiU2, EU2(indu2)+offset,'s','markersize',12,'markerfacecolor',cO,'markeredgecolor','none')
xlim([-psiU0,psiU0])
xticks([-psiU0,psiU0])
xticklabels({'-\Delta\theta/2','\Delta\theta/2'})

x = 5;

end


function plotBumpProfiles(bumpProfiles,optJ,simVars,plotVars)

A   = simVars.A;
th  = simVars.th;
dth = simVars.dth;
ntraj = simVars.ntraj;

cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;

cmap = gray(ntraj+2);

nact = sum(poslin(bumpProfiles));
plot([th;th],[0,A],'--k')
for i = 1:size(bumpProfiles,2)

    activity = bumpProfiles(:,i);
    ii = find(activity>0);

    if optJ
        if i==1
            s = 10;
        else
            s = 6;
        end
        plot(th, activity ,'-','linewidth',lw,'color',cmap(i,:))
        plot(th(ii), activity(ii),'o','markerfacecolor',cmap(i,:),'markeredgecolor','none','markersize',s,'linewidth',lw)
    else
        if nact(i)<mean(nact)
            c = cT;
        else
            c = cO;
        end
        plot(th, activity ,'-','linewidth',lw,'color',c)
        plot(th(ii), activity(ii),'o','markerfacecolor',c,'markeredgecolor','none','markersize',10,'linewidth',lw)
    end
    ylim([0,A])
    xlim([0-dth/2,2*pi-dth/2])
    axis off
end

end

function [actNeurons,indsAct,bumpProfiles,indsBump] = plotMultDriftTraj(JE,N,tmax,tmax2,psi0,optJ,regimes,simVars,plotVars)

% optJ:    binary variable that specifies whether input value of JE is optimal or not
% regimes: binary variable that specifies whether to plot stable/unstable regimes

cK = plotVars.cK;
lw = plotVars.lw;
fs = plotVars.fs;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
ntraj   = simVars.ntraj;
psiSamp = simVars.psiSamp;


jj = find(psi0>=th(3) & psi0<=th(4));
if optJ
    inds = jj(1:2:2*ntraj);
else
    inds = jj([1,ntraj]);
end
indsBump = inds;
bumpProfiles = zeros(N,numel(inds));
m = 1;

vin = 0;
t  = 0:dt:tmax; 
tV = 0:dt:tmax2; 

if regimes
    ii = find(psi0>=2*pi/N & psi0<=pi+2*pi/N);
else
    ii = 1:numel(psi0);
end
actNeurons = zeros(N,numel(ii));
indsAct = ii;


%plot drift trajectories
for i=1:numel(ii)
    [psi,h] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,psi0(ii(i)));
    if i==1
        psi = mod(psi+pi,2*pi)-pi;
    else
        psi = mod(psi,2*pi);
    end
    
    plot(psi,t,'linewidth',lw,'color',cK);    


    xlim([0,2*pi])
    xticks([0,pi,2*pi])
    xticklabels({'0','\pi','2\pi'});   
    
    yticks(0:1:tmax);
    title(['1/J_E = 1/', num2str(round(JE,2))])
    set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');

    actNeurons(poslin(h(1,:))>0,i) = 1;
    if ismember(i,inds)
        bumpProfiles(:,m) = poslin(h(1,:));
        m=m+1;
    end
    
end
        

            
end


function [actNeurons,indsAct,bumpProfiles,indsBump] = plotMultDriftTraj_regimes(JE,JEmin,JEmax,N,tmax,tmax2,psi0,optJ,regimes,simVars,plotVars)

% optJ:    binary variable that specifies whether input value of JE is optimal or not
% regimes: binary variable that specifies whether to plot stable/unstable regimes

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
alphaR = plotVars.alphaR;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
ntraj   = simVars.ntraj;
psiSamp = simVars.psiSamp;


jj = find(psi0>=th(3) & psi0<=th(4));
if optJ
    inds = jj(1:2:2*ntraj);
else
    inds = jj([1,ntraj]);
end
indsBump = inds;
bumpProfiles = zeros(N,numel(inds));
m = 1;

vin = 0;
t  = 0:dt:tmax; 
tV = 0:dt:tmax2; 

if regimes
    ii = find(psi0>=2*pi/N & psi0<=pi+2*pi/N);
else
    ii = 1:numel(psi0);
end
actNeurons = zeros(N,numel(ii));
indsAct = ii;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[dpsiU,dpsiS]   = computeDelNact(tauU,tauS,N);

%plot stable/unstable regimes if nonoptimal J
if ~optJ && regimes

    [cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

    fill([2*pi/N,2*pi/N+dpsiS,2*pi/N+dpsiS,2*pi/N,2*pi/N],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR,'edgecolor','none')
    
    for n=1:N/2
        p0 = -dpsiS+(n+1)*2*pi/N;
        p1 = min(pi+2*pi/N,-dpsiS+(n+1)*2*pi/N+2*dpsiS);
        fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR,'edgecolor','none')

        p0 = dpsiS+n*2*pi/N;
        p1 = min(pi+2*pi/N,dpsiS+n*2*pi/N+2*dpsiU);
        fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR,'edgecolor','none')
    end
end


%plot drift trajectories
for i=1:numel(ii)
    [psi,h] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,psi0(ii(i)));
    if i==1
        psi = mod(psi+pi,2*pi)-pi;
    else
        psi = mod(psi,2*pi);
    end
    istable = find(psi>=2*pi/N & psi<2*pi/N+dpsiS);
    iunstable = [];
    for n=-1:N+1
        p0 = -dpsiS+(n+1)*2*pi/N;
        p1 = min(pi+2*pi/N,-dpsiS+(n+1)*2*pi/N+2*dpsiS);
        istable = [istable;find(psi>=p0 & psi<p1)];
    end
    iunstable = 1:numel(psi);
    iunstable(istable) = [];
        
    %plot drift trajectories
    if i==1
        %psi = mod(psi+pi,2*pi)-pi;
        plot(psi(istable),t(istable),'linewidth',lw,'color',cT);
        plot(psi(istable)+2*pi,t(istable),'linewidth',lw,'color',cT);
        plot(psi(iunstable),t(iunstable),'linewidth',lw,'color',cO);
        plot(psi(iunstable)+2*pi,t(iunstable),'linewidth',lw,'color',cO);
    else
        %psi = mod(psi,2*pi);
        plot(psi(istable),t(istable),'linewidth',lw,'color',cT);
        plot(psi(iunstable),t(iunstable),'linewidth',lw,'color',cO);
    end
    
    if regimes
        xlim([2*pi/N,pi+2*pi/N])
        xticks([2*pi/N,2*pi/N+pi/2,2*pi/N+pi])
        xticklabels({'0','\pi/2','\pi'});   
    else
        xlim([0,2*pi])
        xticks([0,pi,2*pi])
        xticklabels({'0','\pi','2\pi'});   
    end
    
    yticks(0:1:tmax);
    title(['1/J_E = 1/', num2str(round(JE,2))])
    set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');

    actNeurons(poslin(h(1,:))>0,i) = 1;
    if ismember(i,inds)
        bumpProfiles(:,m) = poslin(h(1,:));
        m=m+1;
    end
    
end
        

            
end


function plotSingleDriftTraj(JE,JEmin,JEmax,N,tmax,tmax2,simVars,plotVars)

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
alphaR = plotVars.alphaR;
f = .7;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
dpsi    = simVars.dpsi;
psiSamp = simVars.psiSamp;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,dpsiS]   = computeDelNact(tauU,tauS,N);

[cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

t = 0:dt:tmax;
tV = 0:dt:tmax2;
fill([0,dpsiS,dpsiS,0,0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR)
fill([dpsiS,pi/N,pi/N,dpsiS,dpsiS],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR)

vin = 0;
[psi,~] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,pi/N-dpsi);
psi = mod(psi,2*pi);

i1 = find(psi<dpsiS,1,'first');
i2 = find(psi<dpsi, 1,'first');

plot(psi(1:i1),t(1:i1),'linewidth',lw,'color',cO);
plot(psi(i1:end),t(i1:end),'linewidth',lw,'color',cT);

%plot analytic trajectory
[psiUtmp,tU,psiStmp,tS] = computePsiUnstable(pi/N-dpsi,t,vin,JE,JEmin,JEmax,tau,N);
plot(psiStmp,tS,'--','color',f*cT,'linewidth',lw);
plot(psiUtmp,tU,'--','color',f*cO,'linewidth',lw);

plot([0,pi/N],[t(i2),t(i2)],'--k');
plot([dpsi,dpsi],[0,tmax],'--k');
plot([pi/N-dpsi,pi/N-dpsi],[0,tmax],'--k');

xlim([0,pi/N])
xticks([0,pi/N])
xticklabels({'0','\pi/N'});    

set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');
end

function plotDoubleDrivenTraj(vin,JE,JEmin,JEmax,N,tmax,tmax2,simVars,plotVars)

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
ms = plotVars.ms;
alphaR = plotVars.alphaR;
f = .7;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

[lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,dpsiS] = computeDelNact(tauU,tauS,N);
psiU = computeUnstableFPs(vin,lambdaU,N); 
psiS = computeStableFPs(  vin,lambdaS); 
vmin = computeVmin(tauS,dpsiS);

[cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

t  = 0:dt:tmax;
tV = 0:dt:tmax2;
fill([0,dpsiS,dpsiS,0,0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR)
fill([dpsiS,pi/N,pi/N,dpsiS,dpsiS],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR)


[psi,~] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,0);
psi = mod(psi,2*pi);
if vin>vmin
    i1 = find(psi>dpsiS,1,'first');
    plot(psi(1:i1-1),t(1:i1-1),'linewidth',lw,'color',cT);
    plot(psi(i1:end),t(i1:end),'linewidth',lw,'color',cO);
else
    psi = mod(psi+pi,2*pi)-pi;
    plot(psi,t,'linewidth',lw,'color',cT);
    
    %plot analytic trajectory
    psiStmp = computePsiStable(t,vin,JE,JEmin,JEmax,tau,N);
    plot(psiStmp,t,'--','color', f*cT,'linewidth',lw);
end

if vin<=vmin
    [psi,~] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,psiU-(psiU-dpsiS)./10);
    psi = mod(psi,2*pi);
    i1 = find(psi<dpsiS,1,'first');
    plot(psi(1:i1),t(1:i1),'linewidth',lw,'color',cO);
    plot(psi(i1:end),t(i1:end),'linewidth',lw,'color',cT);
    
    %plot analytic trajectory
    psiStmp = computePsiStable(t,vin,JE,JEmin,JEmax,tau,N);
    plot(psiStmp,t,'--','color', f*cT,'linewidth',lw);
    
    [psiUtmp,tU,psiStmp,tS] = computePsiUnstable(psiU-(psiU-dpsiS)./10,t,vin,JE,JEmin,JEmax,tau,N);
    plot(psiStmp,tS,'--','color',f*cT,'linewidth',lw);
    plot(psiUtmp,tU,'--','color',f*cO,'linewidth',lw);
        
else
    %plot analytic trajectory
    [~,~,tS,~, psiStail,psiUtail,tStail,tUtail] = computePsi(vin,JE,JEmin,JEmax,tau,N);
    plot(psiStail,tStail,'--','linewidth',lw,'color',f*cT);
    plot(psiUtail,tUtail+tS(end),'--','linewidth',lw,'color',f*cO);
end

%plot FPs
plot([psiU,psiU],[0,tmax],'--k');
plot([psiS,psiS],[0,tmax],'--k');

plot(psiS,tmax(end),'s','markersize',ms,'markerfacecolor','k','markeredgecolor','none')
plot(psiU,0        ,'s','markersize',ms,'markerfacecolor','w','markeredgecolor','k')

ylim([0,tmax(end)])
xlim([0,pi/N])
xticks([0,pi/N])
xticklabels({'0','\pi/N'});    
ylabel('time (s)')
xlabel('bump orientation')

set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');
end

function plotHeatmap(JEmin,JEmax,N,simVars,plotVars)

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;

tau = simVars.tau;


nsJ = 200;
nsW = 1000;
psiVec = linspace(0,pi/N,nsW);

JEint = 1./linspace(1./JEmax,1./JEmin,nsJ);

yt = linspace(1./JEmax,1./JEmin,7);
for i=1:7
    yl{i} = ['1/',num2str(round(1/yt(i),2))];
end

%compute drift rates
[lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
[~,dpsiS] = computeDelNact(tauU,tauS,N);

%build matrix
M = [];
for i=1:nsJ
    ii = find(psiVec>dpsiS(i),1,'first');
    M = [M;[abs(lambdaS(i))*ones(1,numel(1:ii-1)),lambdaU(i)*ones(1,numel(ii:nsW))]];
end

%build colormap
nC = 200;
f  = max(abs(lambdaS))./max(abs(lambdaU));
cmapU = buildColormap(cO,cW,nC);
cmapS = buildColormap(cT,cW,floor(f*nC));
npad  = nC-floor(f*nC);
cmapS = [cmapS;repmat(cmapS(end,:),[npad,1])];

cmap = [flipud(cmapS);cmapU(2:end,:)];

%plot results
imagesc(psiVec,1./JEint,M);colormap(flipud(gray));caxis([0,max(lambdaU)])
hold on;
plot(dpsiS,1./JEint,'-k','linewidth',lw)
set(gca,'ydir','normal','fontsize',fs,'TickDir','out');
xticks([0,pi/N])
xticklabels({'0','\pi/N'})
yticks(linspace(1./JEmax,1./JEmin,7))
yticklabels(yl)
xlabel('bump orientation')
ylabel('excitation 1/J_E')
    
end

function plotHeatmap_color(JEmin,JEmax,N,simVars,plotVars)

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;

tau = simVars.tau;


nsJ = 200;
nsW = 1000;
psiVec = linspace(0,pi/N,nsW);

JEint = 1./linspace(1./JEmax,1./JEmin,nsJ);

yt = linspace(1./JEmax,1./JEmin,7);
for i=1:7
    yl{i} = ['1/',num2str(round(1/yt(i),2))];
end

%compute drift rates
[lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JEint,JEmin,JEmax,tau);
[~,dpsiS] = computeDelNact(tauU,tauS,N);

%build matrix
M = [];
for i=1:nsJ
    ii = find(psiVec>dpsiS(i),1,'first');
    M = [M;[(lambdaS(i))*ones(1,numel(1:ii-1)),lambdaU(i)*ones(1,numel(ii:nsW))]];
end
M = M';
%build colormap
nC = 200;
f  = max(abs(lambdaS))./max(abs(lambdaU));
cmapU = buildColormap(cO,cW,nC);
cmapS = buildColormap(cT,cW,floor(f*nC));
npad  = nC-floor(f*nC);
cmapS = [cmapS;repmat(cmapS(end,:),[npad,1])];

cmap = [flipud(cmapS);cmapU(2:end,:)];

%plot results
imagesc(1./JEint,psiVec,M);colormap(flipud(cmap));caxis([-max(lambdaU),max(lambdaU)])
hold on;
plot(1./JEint,dpsiS,'-k','linewidth',lw)
set(gca,'ydir','normal','fontsize',fs,'TickDir','out');
yticks([0,pi/N])
yticklabels({'0','\pi/N'})
xticks(linspace(1./JEmax,1./JEmin,7))
xticklabels(yl)
ylabel('bump orientation')
xlabel('excitation 1/J_E')
    
end

function plotMultDrivenTrajBelowThresh(vin,JE,JEmin,JEmax,N,tmax,tmax2,simVars,plotVars,psiInit)

if nargin<9
    psiInit = 0;
end

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
alphaR = plotVars.alphaR;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,dpsiS]   = computeDelNact(tauU,tauS,N);
    
[cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

t  = 0:dt:tmax;
tV = 0:dt:tmax2;
fill([0,dpsiS,dpsiS,0,0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR)
fill([dpsiS,pi/N,pi/N,dpsiS,dpsiS],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR)

cmap = gray(numel(vin)+2);
cmap = flipud(cmap(1:end-2,:));
for v = 1:numel(vin)
    [psi,~] = simPsi(JE,vin(v),th,A,C0,psiSamp,tau,N,t,tV,del,psiInit);
    psi = unwrap(psi);
    plot(psi,t,'linewidth',lw,'color',cmap(v,:))
end


xlim([0,pi/N])
xticks([0,pi/N])
xticklabels({'0','\pi/N'});    

set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');
end

function plotMultDrivenTrajAboveThresh_regimes(vin,JE,JEmin,JEmax,N,tmax,tmax2,optJ,regimes,simVars,plotVars,psiInit)

% optJ:    binary variable that specifies whether input value of JE is optimal or not
% regimes: binary variable that specifies whether to plot stable/unstable regimes


if nargin<12
    psiInit = 0;
end

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
alphaR = plotVars.alphaR;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

t  = 0:dt:tmax;  
tV = 0:dt:tmax2;  

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[dpsiU,dpsiS]   = computeDelNact(tauU,tauS,N);

%plot stable/unstable regimes if nonoptimal J
if ~optJ && regimes

    [cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

    fill([0,dpsiS,dpsiS,0,0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR,'edgecolor','none')
    for n=1:N/2
        p0 = -dpsiS+n*2*pi/N;
        p1 = min(pi,-dpsiS+n*2*pi/N+2*dpsiS);
        fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR,'edgecolor','none')

        p0 = dpsiS+(n-1)*2*pi/N;
        p1 = min(pi,dpsiS+(n-1)*2*pi/N+2*dpsiU);
        fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR,'edgecolor','none')
    end
end

cmap = gray(numel(vin)+2);
cmap = flipud(cmap(1:end-2,:));

alphas = linspace(0,1,numel(vin)+2);
alphas = alphas(3:end);


%plot driven trajectories
for i=1:numel(vin)
    [psi,~] = simPsi(JE,vin(i),th,A,C0,psiSamp,tau,N,t,tV,del,psiInit);
    psi = mod(psi,2*pi);
    
    istable = find(psi>=2*pi/N & psi<2*pi/N+dpsiS);
    iunstable = [];
    for n=-N:2*N+1
        p0 = -dpsiS+(n+1)*2*pi/N;
        p1 = min(pi+2*pi/N,-dpsiS+(n+1)*2*pi/N+2*dpsiS);
        istable = [istable;find(psi>=p0 & psi<p1)];
    end
    iunstable = 1:numel(psi);
    iunstable(istable) = [];
    istable = sort(istable);
    iunstable = sort(iunstable);
    

    ix = find(diff(istable)>1);
    if numel(ix)>0
        plot(psi(istable(1:ix(1))),t(istable(1:ix(1))),'linewidth',lw,'color',[cT,alphas(i)]);
        for j=1:numel(ix)-1
            plot(psi(istable(ix(j)+1:ix(j+1))),t(istable(ix(j)+1:ix(j+1))),'linewidth',lw,'color',[cT,alphas(i)]);
        end
        plot(psi(istable(ix(end)+1:end)),t(istable(ix(end)+1:end)),'linewidth',lw,'color',[cT,alphas(i)]);

    else
        plot(psi(istable),t(istable),'linewidth',lw,'color',[cT,alphas(i)]);
    end
    ix = find(diff(iunstable)>1);
    if numel(ix)>0
        plot(psi(iunstable(1:ix(1))),t(iunstable(1:ix(1))),'linewidth',lw,'color',[cO,alphas(i)]);
        for j=1:numel(ix)-1
            plot(psi(iunstable(ix(j)+1:ix(j+1))),t(iunstable(ix(j)+1:ix(j+1))),'linewidth',lw,'color',[cO,alphas(i)]);
        end
        plot(psi(iunstable(ix(end)+1:end)),t(iunstable(ix(end)+1:end)),'linewidth',lw,'color',[cO,alphas(i)]);
    else
        plot(psi(iunstable),t(iunstable),'linewidth',lw,'color',[cO,alphas(i)]);
    end
    
    if regimes
        xlim([0,pi])
        xticks([0,pi/2,pi])
        xticklabels({'0','\pi/2','\pi'});   
    else
        xlim([0,2*pi])
        xticks([0,pi,2*pi])
        xticklabels({'0','\pi','2\pi'});   
    end
    
    yticks(0:1:tmax);
    title(['1/J_E = 1/', num2str(round(JE,2))])
    set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');

    
end
        

            
end


function plotMultDrivenTrajAboveThresh(vin,JE,N,tmax,tmax2,simVars,plotVars,psiInit)

% optJ:    binary variable that specifies whether input value of JE is optimal or not
% regimes: binary variable that specifies whether to plot stable/unstable regimes

if nargin<8
    psiInit = 0;
end

lw = plotVars.lw;
fs = plotVars.fs;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;


t  = 0:dt:tmax;  
tV = 0:dt:tmax2;  

cmap = gray(numel(vin)+2);
cmap = flipud(cmap(1:end-2,:));

%plot driven trajectories
for i=1:numel(vin)
    [psi,~] = simPsi(JE,vin(i),th,A,C0,psiSamp,tau,N,t,tV,del,psiInit);
    psi = mod(psi,2*pi);

    plot(psi,t,'linewidth',lw,'color',cmap(i,:));
    
    yticks(0:1:tmax);
    title(['1/J_E = 1/', num2str(round(JE,2))])
    set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');

end
       
            
end



function [psiUout,tUout,psiSout,tSout] = plotSingleDrivenTrajAboveThresh(vin,JE,JEmin,JEmax,N,tmax,tmax2,simVars,plotVars,psiInit)

if nargin<9
    psiInit = 0;
end

cW = plotVars.cW;
cO = plotVars.cO;
cT = plotVars.cT;
lw = plotVars.lw;
fs = plotVars.fs;
alphaR = plotVars.alphaR;

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[dpsiU,dpsiS]   = computeDelNact(tauU,tauS,N);
    
[cS,cU] = getDriftColor(cW,cT,cO,JE,JEmin,JEmax,tau,1000);

t  = 0:dt/5:tmax;
tV = 0:dt:tmax2;

rangeS = [0,dpsiS];
rangeU = [dpsiS,dpsiS+2*dpsiU];
fill([0,dpsiS,dpsiS,0,0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR)
fill([dpsiS,dpsiS+2*dpsiU,dpsiS+2*dpsiU,dpsiS,dpsiS],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR)
for n=1
    p0 = -dpsiS+n*2*pi/N;
    p1 = min(pi/2,-dpsiS+n*2*pi/N+2*dpsiS);
    fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cS,'FaceAlpha',alphaR)
    rangeS = [rangeS;[p0,p1]];
    
    p0 = dpsiS+n*2*pi/N;
    p1 = min(pi/2,dpsiS+n*2*pi/N+2*dpsiU);
    fill([p0,p1,p1,p0,p0],[0,0,tmax,tmax,0],cU,'FaceAlpha',alphaR)
    rangeU = [rangeU;[p0,p1]];
end

[psi,~] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,tV,del,psiInit);
psi = unwrap(psi);

[vmax,vmin] = computeBumpVelocity(psi,t,dpsiS,N);
for i=1:size(rangeS,1)
    indsS = find(psi>=rangeS(i,1) & psi<=rangeS(i,2));
    plot(psi(indsS),t(indsS),'linewidth',lw,'color',cT)
    
    if numel(indsS)>0
        %plot min velocity as bump leaves the stable regime
        i0 = indsS(end);
        b = t(i0)-psi(i0)./vmin;
        psiS = linspace(psi(i0)-(pi/N).*vmin, psi(i0)+(pi/N).*vmin);
        tS = psiS./vmin+b;
        plot(psiS,tS,'-k','linewidth',lw)
        
        if i==1
            psiSout = psiS-psi(i0);
            tSout = tS-t(i0);
        end
    end
end
for i=1:size(rangeU,1)
    indsU = find(psi>=rangeU(i,1) & psi<=rangeU(i,2));
    plot(psi(indsU),t(indsU),'linewidth',lw,'color',cO)
    
    if numel(indsU)>0
        %plot max velocity as bump leaves the unstable regime
        i0 = indsU(end);
        b = t(i0)-psi(i0)./vmax;
        psiU = linspace(psi(i0)-(pi/N).*vmax/3, psi(i0)+(pi/N).*vmax/3);
        tU = psiU./vmax+b;
        plot(psiU,tU,'-k','linewidth',lw)
        
        if i==1
            psiUout = psiU-psi(i0);
            tUout = tU-t(i0);
        end
    end
end



ylim([0,tmax/2])
xlim([0,pi/2])
xticks([0,pi/4,pi/2])
xticklabels({'0','\pi/4','\pi/2'});    

set(gca,'ydir','reverse','fontsize',fs,'TickDir','out');
end

function [M,Jint,maxM,dn] = plotNdepHeatmap(var,nSamples,simVars,plotVars,fighandle,subplotinds)

fs = plotVars.fs;

tau = simVars.tau;
dpsi = simVars.dpsi;

if strcmp(var,'drift')
    %cax = [0,.5];
    cax = [-5,0];
    fun = @max;
elseif strcmp(var,'vmin')
    %cax = [0,1];
    cax = [-5,0];
    fun = @max;
elseif strcmp(var,'lin')
    %cax = [.5,1];
    cax = [-5,0];
    fun = @min;
else
    error('unrecognized input variable')
end

%linearly spaced in 1/JE
nS   = 100000;
Jint = linspace(0,.5,nS);

M    = nan(numel(Jint)  ,numel(nSamples));
maxM = nan(max(nSamples),numel(nSamples));

for n = 1:numel(nSamples)
    [JEopt,~,Nact] = computeOptJE(nSamples(n));
    JEopt = fliplr(JEopt);
    JEopt = [10000,JEopt,2];
    
    Nact  = fliplr(Nact);
    Nact  = [1,Nact,nSamples(n)-1];
        
    for j=2:numel(JEopt)-2%j=1:numel(JEopt)-1
        inds  = find(Jint>1./JEopt(j) & Jint<=1./JEopt(j+1));
        JEinv = Jint(inds);
        JE = 1./JEinv;
        
        
        [~,~,tauU, tauS ] = computeTimeConstants(JE, JEopt(j+1),JEopt(j),tau);
        if strcmp(var,'drift')
            tDrift     = computeDriftTime(tauU, tauS);
            driftRate  = computeNetDriftRate(tDrift, nSamples(n));
            vec  = driftRate;   
        else
            [~,dpsiS]  = computeDelNact(tauU, tauS, nSamples(n));
            vmin  = computeVmin(tauS, dpsiS);
            if strcmp(var,'vmin')
                vec  = vmin;
            elseif strcmp(var,'lin')
                vec  = (1.5-vmin )./(1.5+vmin);
            end
        end    
        
        dn = Nact(j)-nSamples(n)/2 + max(nSamples)/2+1;
        M( inds,n) = vec;
        maxM(dn,n) = fun(vec);

            
    end
end
subplot(subplotinds(1),subplotinds(2),subplotinds(3),'Parent',fighandle);
h = imagesc(nSamples,Jint,log(M),'alphadata',~isnan(M));caxis(cax)
hold on;
set(gca,'fontsize',fs,'ydir','normal','tickdir','out')
ylim([0,.5])

for n = 1:numel(nSamples)
    [JEopt,~,~] = computeOptJE(nSamples(n));
    plot(nSamples(n),1./JEopt,'o','markerfacecolor','none','markeredgecolor','r','markersize',8)
end
ylabel('local excitation 1/J_E');
xticks(nSamples)
xlabel('size N')
    
end


function [M,Jint,maxM,dn] = plotNdepHeatmapScratch(var,nSamples,simVars,plotVars,fighandle,subplotinds)

fs = plotVars.fs;

tau = simVars.tau;
dpsi = simVars.dpsi;

if strcmp(var,'drift')
    %cax = [0,.5];
    cax = [-5,0];
    fun = @max;
elseif strcmp(var,'vmin')
    %cax = [0,1];
    cax = [-5,0];
    fun = @max;
elseif strcmp(var,'lin')
    %cax = [.5,1];
    cax = [-5,0];
    fun = @min;
else
    error('unrecognized input variable')
end

%linearly spaced in 1/JE
nS   = 100000;
Jint = linspace(0,.5,nS);

M    = nan(numel(Jint)  ,numel(nSamples));
maxM = nan(max(nSamples),numel(nSamples));

for n = 1:numel(nSamples)
    [JEopt,~,Nact] = computeOptJE(nSamples(n));
    JEopt = fliplr(JEopt);
    JEopt = [10000,JEopt,2];
    
    Nact  = fliplr(Nact);
    Nact  = [1,Nact,nSamples(n)-1];
        
    for j=1:numel(JEopt)-1
        inds  = find(Jint>1./JEopt(j) & Jint<=1./JEopt(j+1));
        JEinv = Jint(inds);
        JE = 1./JEinv;
        
        
        [~,~,tauU, tauS ] = computeTimeConstants(JE, JEopt(j+1),JEopt(j),tau);
        if strcmp(var,'drift')
            tDrift     = computeDriftTime(tauU, tauS);
            driftRate  = computeNetDriftRate(tDrift, nSamples(n));
            vec  = driftRate;   
        else
            [~,dpsiS]  = computeDelNact(tauU, tauS, nSamples(n));
            vmin  = computeVmin(tauS, dpsiS);
            if strcmp(var,'vmin')
                vec  = vmin;
            elseif strcmp(var,'lin')
                vec  = (1.5-vmin )./(1.5+vmin);
            end
        end    
        
        dn = Nact(j)-nSamples(n)/2 + max(nSamples)/2+1;
        M( inds,n) = vec;
        maxM(dn,n) = fun(vec);

            
    end
end
subplot(subplotinds(1),subplotinds(2),subplotinds(3),'Parent',fighandle);
h = imagesc(nSamples,Jint,log(M),'alphadata',~isnan(M));caxis(cax);colormap(viridis)
hold on;
set(gca,'fontsize',fs,'ydir','normal','tickdir','out')
ylim([0,.5])

for n = 1:numel(nSamples)
    [JEopt,~,~] = computeOptJE(nSamples(n));
    plot(nSamples(n),1./JEopt,'o','markerfacecolor','none','markeredgecolor','r','markersize',8)
end
ylabel('local excitation 1/J_E');
xticks(nSamples)
xlabel('size N')
    
end

function [dJ,J] = plotNdepThresh(var,thresh,nSamples,simVars,plotVars,fighandle,subplotinds)

lw = plotVars.lw;
fs = plotVars.fs;
cR = plotVars.cR;
cmap = plotVars.cmapN;

tau = simVars.tau;

if strcmp(var,'drift')
    c = (exp(1)-1)./(2*exp(1));
    f = tau./(c.*pi);
    xlm = [-7,-1];
elseif strcmp(var,'vmin')
    f = 2.*tau./pi;
    xlm = [-7,-1];
elseif strcmp(var,'lin')
    f = 1.5*tau./pi; %1.5 is the input velocity used above
    xlm = [-10,-4];
else
    error('unrecognized input variable')
end

nS = numel(nSamples);
dJ = nan(max(nSamples)-2,nS);
J  = nan(max(nSamples)-2,nS);


for n = 1:nS

    [JEopt,~,Nact] = computeOptJE(nSamples(n));   
    Jr(n) = max(JEopt)-min(JEopt);
    for j=1:numel(JEopt)
        
        if j>1 && j<numel(JEopt)
            J1  = 1./((1./JEopt(j)+1./JEopt(j-1))/2);
            J2  = 1./((1./JEopt(j+1)+1./JEopt(j))/2);
            JE1 = JEopt(j-1):.000001:JEopt(j);JE1 = JE1(1:end-1);
            JE2 = JEopt(j):.000001:JEopt(j+1);
            JE  = [JE1,JE2];

            [~,~,tauU1,tauS1] = computeTimeConstants(JE1,JEopt(j-1),JEopt(j),  tau);
            [~,~,tauU2,tauS2] = computeTimeConstants(JE2,JEopt(j),  JEopt(j+1),tau);
            tauU = [tauU1,tauU2];
            tauS = [tauS1,tauS2];
        elseif j==1
            J1  = JEopt(j);
            J2  = 1./((1./JEopt(j+1)+1./JEopt(j))/2);
            JE = JEopt(j):.000001:JEopt(j+1);
            [~,~,tauU,tauS] = computeTimeConstants(JE,JEopt(j),JEopt(j+1),tau);
        elseif j==numel(JEopt)
            J1  = 1./((1./JEopt(j)+1./JEopt(j-1))/2);
            J2  = JEopt(j);
            JE = JEopt(j-1):.000001:JEopt(j);
            [~,~,tauU,tauS] = computeTimeConstants(JE,JEopt(j-1),JEopt(j),tau);
        end
            

        if strcmp(var,'drift')
            tDrift  = computeDriftTime(tauU, tauS);
            vec  = computeNetDriftRate(tDrift, nSamples(n)); 
        else
            [~,dpsiS] = computeDelNact(tauU,tauS,nSamples(n));
            vmin = computeVmin(tauS,dpsiS);
            if strcmp(var,'vmin')
                vec  = vmin;
            elseif strcmp(var,'lin')
                vec  = (1.5-vmin)./(1.5+vmin);
            end
        end
        
        [~,i1] = min(abs(JE-J1));
        [~,i2] = min(abs(JE-J2));
        JEint = JE( i1:i2);
        vint  = vec(i1:i2);
        dn = Nact(j)-nSamples(n)/2 + max(nSamples)/2+1;
        
        if strcmp(var,'lin')
            j1 = find(vint>thresh,1,'first');
            j2 = find(vint>thresh,1,'last' );
        else
            j1 = find(vint<thresh,1,'first');
            j2 = find(vint<thresh,1,'last' );
        end

        if j>1 && j<numel(JEopt)
            dJ(dn,n) = (JEint(j2)-JEint(j1));
        else
            dJ(dn,n) = 2*(JEint(j2)-JEint(j1));
        end
        J( dn,n) = JEopt(j);
        
        
        if n==1 && j==2
            ymax = [2,.3];
            xmin = [J1,JEopt(j)-.2];
            xmax = [J2,JEopt(j)+.2];
            for k=1:2
                subplot(subplotinds(1),subplotinds(2),subplotinds(3)+k-1,'Parent',fighandle);hold on;
                inds1 = 1:100:numel(JE1);
                inds2 = 1:100:numel(JE2);

                fscale = 100;
                jj1 = find(vint<fscale*thresh,1,'first');
                jj2 = find(vint<fscale*thresh,1,'last' );

                m = 2./(nSamples(n).*JEopt(j).*f);
                plot(JE,vec,'-k','linewidth',lw)
                plot(JE1(inds1),-m*(JE1(inds1)-JEopt(j)))
                plot(JE2(inds2),m*(JE2(inds2)-JEopt(j)))
                plot([J1,J2],fscale*[thresh,thresh],'--k')
                plot([JEint(jj1),JEint(jj2)],[0,0],'-k','linewidth',4)
                plot([JEint(jj1),JEint(jj2)],[0,0],'-k','linewidth',4)
                plot([JEopt(j),JEopt(j)],[0,ymax(k)],'--','color',cR)
                plot([JEint(jj1),JEint(jj1)],[0,fscale*thresh],'--k')
                plot([JEint(jj2),JEint(jj2)],[0,fscale*thresh],'--k')

                xlim([xmin(k),xmax(k)])
                ylim([0,ymax(k)]);
                set(gca,'fontsize',fs,'tickdir','out')
                xlabel('J_E')
            end
            
        end
        
    end
end

subplot(subplotinds(1),subplotinds(2),subplotinds(3)+2,'Parent',fighandle);hold on;
jj = computeOptJE(max(nSamples));
je = 1./linspace(1./max(jj),1./min(jj),10000);
for i=1:numel(nSamples)
    tol = nSamples(i).*thresh.*je.*f;
    plot(je,tol,'-','color',.9*cmap(i,:))
    plot(J(:,i),dJ(:,i),'ok','markerfacecolor',cmap(i,:),'markersize',12)
end
set(gca,'fontsize',fs,'tickdir','out','YScale', 'log','XScale', 'log')



[JEall,~,~] = computeOptJE(4);
subplot(subplotinds(1),subplotinds(2),subplotinds(3)+3,'Parent',fighandle);hold on;
tol = nSamples.*thresh.*JEall.*f;
plot(nSamples,tol,'-k','linewidth',lw)
for i=1:numel(nSamples)
    plot(nSamples(i),dJ(11,i),'ok','markerfacecolor',cmap(i,:),'markersize',12)
end
set(gca,'fontsize',fs,'tickdir','out','YScale', 'log','XScale', 'log')
xticks(nSamples)
xlabel('size N')
ylabel('tolerance')

subplot(subplotinds(1),subplotinds(2),subplotinds(3)+4,'Parent',fighandle);hold on;
nn = repmat(nSamples,[size(J,1),1]);
nn(isnan(J)) = nan;
tol = nn.*J.*f.*thresh;
tolLB = nSamples.^2./(1-cos(2.*pi./nSamples)).*f.*thresh;
plot(nSamples,tolLB,'--r','linewidth',lw)
plot(nSamples,nansum(tol  ),'-ok','linewidth',lw)
plot(nSamples,nansum(dJ   ),'-k','linewidth',lw)
set(gca,'fontsize',fs,'tickdir','out')
xticks(nSamples)
xlabel('size N')
ylabel('total tolerance')


    
end

function errVar = computeNoiseVar(N,JE,t,psi0,th,nSims,noiseLevel,A,C0,tau,del,psiSamp)

err = [];
vin = 0;

parfor i=1:nSims
    aNoise = normrnd(0,noiseLevel,N,length(t));
    [psi,~] = simPsi(JE,vin,th,A,C0,psiSamp,tau,N,t,[],del,psi0,aNoise);

    %unwrap trajectory
    dpsi = diff(psi);
    tshift = find(abs(dpsi)>pi);
    psiU = psi;
    for j=1:numel(tshift)
        psiU(tshift(j)+1:end) = psiU(tshift(j)+1:end)-sign(dpsi(tshift(j)))*2*pi;
    end

    err = [err;(psiU'-psi0).^2];

end
errVar = sum(err)./(nSims-1);
end

function genNoiseVar_varJE(N,noiseLevel,simVars)

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
psiSamp = simVars.psiSamp;


%generate noise tolerance for different JE
nSims = 10000;
tMax = 20;
t = 0:dt:tMax;
tS = find(t>10,1,'first');

[JEopt,~,~] = computeOptJE(N);
th  = linspace(0, 2*pi, N+1); 
th(end) = [];
psi0 = 0;

filename = ['noiseSimulations/noise_varJE_N',num2str(N),'_noise',num2str(round(noiseLevel,3)),'.mat'];
errVar = [];
for k = 1:numel(JEopt)
    JE = JEopt(k);
    eV = computeNoiseVar(N,JE,t,psi0,th,nSims,noiseLevel,A,C0,tau,del,psiSamp);
    p = polyfit(t(tS:end),eV(tS:end),1);
    errTol(k) = p(1);
    errVar = [errVar;eV];
    save(filename,'JEopt','errVar','errTol','nSims','t','noiseLevel')
end

end

function plotNoiseVar_varJE(N,noiseLevel,plotVars)
lw = plotVars.lw;
fs = plotVars.fs;
cmap = plotVars.cmapN;


filename = ['noiseSimulations/noise_varJE_N',num2str(N),'_noise',num2str(round(noiseLevel,3)),'.mat'];
load(filename)

figure;
set(gcf,'Position',[200 200 800 800],'color','w');

subplot(2,2,1);
plot(t,errVar,'linewidth',lw);
legend({['JE = ',num2str(JEopt(1))],['JE = ',num2str(JEopt(2))],['JE = ',num2str(JEopt(3))]})
set(gca,'fontsize',fs)
xlabel('time')
ylabel('error variance')

subplot(2,2,3);hold on;
tS = find(t>10,1,'first');
plot(t(tS:end),errVar(:,tS:end))
for i=1:3
    p = polyfit(t(tS:end),errVar(i,tS:end),1);
    plot(t(tS:end),p(1).*t(tS:end)+p(2),'--k')
end
set(gca,'fontsize',fs,'YScale', 'log','XScale', 'log')
xlabel('time')
ylabel('error variance')

subplot(2,2,2);hold on
N = [6,8,10,12,14];
for i=1:numel(N)
    filename = ['noiseSimulations/noise_varJE_N',num2str(N(i)),'_noise',num2str(round(noiseLevel,3)),'.mat'];
    load(filename)
    plot(JEopt,1./errTol,'-o','color',cmap(i+1,:))
end
set(gca,'fontsize',fs,'YScale', 'log','XScale', 'log')
xlabel('optimal JE')
ylabel('robustness')

end

function plotNoiseVar_varN(JE,noiseLevels,plotVars)
lw = plotVars.lw;
fs = plotVars.fs;
cmap = plotVars.cmapN;

figure;
set(gcf,'Position',[200 200 800 800],'color','w');hold on
for i=1:numel(noiseLevels)
    filename = ['noiseSimulations/noise_varN_JE',num2str(JE),'_noise',num2str(round(noiseLevels(i),3)),'.mat'];
    load(filename)
    subplot(2,2,1);hold on;
    plot(nN,1./errTol,'-k','linewidth',lw);
    scatter(nN,1./errTol,100,nN,'filled');colormap(cmap)

    subplot(2,2,2);hold on;
    plot(nN,1./errTol,'-k','linewidth',lw);
    scatter(nN,1./errTol,100,nN,'filled');colormap(cmap)

    p = polyfit(nN,1./errTol,1);
    a(i) = p(1);
    b(i) = p(2); 
end
subplot(2,2,1);hold on;
xlim([2,22])
set(gca,'fontsize',fs,'YScale', 'log','XScale', 'log')
xlabel('network size N')
ylabel('noise robustness')

subplot(2,2,2);hold on;
xlim([2,22])
set(gca,'fontsize',fs)
xlabel('network size N')
ylabel('noise robustness')

options = optimoptions('fminunc','Display','off');
x = 1./(noiseLevels.^2);
subplot(2,2,3);hold on;
plot(x,a,'o')
a0 = fminunc(@(m)sum((m.*x-a).^2),0,options);
plot(x,a0.*x,'--k')
set(gca,'fontsize',fs,'YScale', 'log','XScale', 'log')
xlabel('1/\sigma^2')
ylabel('slope')
title(['a = ',num2str(round(a0,3)),'/\sigma^2'])

subplot(2,2,4);hold on;
plot(x,b,'o')
b0 = fminunc(@(m)sum((m.*x-b).^2),0,options);
plot(x,b0.*x,'--k')
set(gca,'fontsize',fs)
xlabel('1/\sigma^2')
ylabel('offset')
title(['b = ',num2str(round(b0,3)),'/\sigma^2'])


for i=1:numel(noiseLevels) 
    filename = ['noiseSimulations/noise_varN_JE',num2str(JE),'_noise',num2str(round(noiseLevels(i),3)),'.mat'];
    load(filename)
    subplot(2,2,1);hold on;
    plot(nN,(a0.*nN+b0)./(noiseLevels(i).^2),'--k')
end
title('(aN+b)/\sigma^2')

end


function genNoiseVar_varN(JE,noiseLevel,simVars)

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
psiSamp = simVars.psiSamp;


%generate noise tolerance for different JE
nSims = 10000;
tMax = 20;
t = 0:dt:tMax;
tS = find(t>10,1,'first');

psi0 = 0;
nN = 4:2:20;
errVar = [];
filename = ['noiseSimulations/noise_varN_JE',num2str(JE),'_noise',num2str(round(noiseLevel,3)),'.mat'];
for k = 1:numel(nN)
    N = nN(k);
    th  = linspace(0, 2*pi, N+1); 
    th(end) = [];
    eV = computeNoiseVar(N,JE,t,psi0,th,nSims,noiseLevel,A,C0,tau,del,psiSamp);
    p = polyfit(t(tS:end),eV(tS:end),1);
    errTol(k) = p(1);
    errVar = [errVar;eV];
    save(filename,'nN','errVar','errTol','nSims','t','noiseLevel')
end




end

%----------------------- ANALYTIC READOUTS -------------------------------%

function [Ws,Wa,JI] = getConnectivity(N,JE,A,C0,psiSamp,del,scale)

if nargin<7
    scale = 1;
end

% get preferred headings
th = linspace(0,2*pi,N+1);th(end) = [];

% get inhibition
if scale
    JI = computeJI(JE,A,C0,psiSamp,N);
else
    JI = -100;
end

% get weight matrices
Ws = 1/N*(JI + JE*cos(th - th'));
Wa = 1/N*cos(th - th' + del);
end

function W = getFullConnectivity(Ws,Wa,v,tau)
W = (Ws+v.*Wa-eye(size(Ws,1)))./tau;
end

function [JEopt,JEinv,Nact] = computeOptJE(N)
Nact = fliplr(2:N-2);
n = Nact-N/2;
for i=1:numel(n)
    JEinv(i) = 1./4 + 1./(2.*N) .*( n(i) + sin(2*pi*n(i)./N)./sin(2*pi./N) );
end
JEopt = 1./JEinv;
end

function [lambdaU,lambdaS,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau)
lambdaU = (JE./JEmin-1)./tau;   % drift rate in unstable regime 
lambdaS = (JE./JEmax-1)./tau;   % drift rate in stable regime
tauU = 1./lambdaU;
tauS = 1./abs(lambdaS);
end

function [dpsiU,dpsiS] = computeDelNact(tauU,tauS,N)
dpsiU = pi./N.*(1./(1+tauS./tauU)); % half-width of stable regime
dpsiS = pi./N.*(1./(1+tauU./tauS)); % half-width of unstable regime
end

function tdrift = computeDriftTime(tauU,tauS)
% time to drift from psi = pi/N-epsU to psi = epsS;
tdrift = tauU+tauS;
end

function lambdaD = computeNetDriftRate(tdrift,N)
dTheta = 2*pi/N;
c = (exp(1)-1)./(2*exp(1));
lambdaD = c.*dTheta./tdrift;
end

function vmin = computeVmin(tauS,dpsiS)
% minimum integrable velocity
vmin = dpsiS./(tauS);
end

function r = computeLinearity(v,vmin)
% linearity of integration
r = (v-vmin)./(v+vmin);
end

function psiS = computeStableFPs(v,lambdaS)
psiS = v./abs(lambdaS);
end

function psiU = computeUnstableFPs(v,lambdaU,N)
psiU = pi/N-v./abs(lambdaU);
end

function [evU,evS] = computeEvals(JE,JEmin,JEmax,tau)
evU = (JE./JEmin-1)./tau;
evS = (JE./JEmax-1)./tau;
end


%----------------------- SIMULATION READOUTS -----------------------------%

function vmin = simVmin(JE,JEmin,JEmax,N,tmax2,simVars)

A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[dpsiU,dpsiS] = computeDelNact(tauU,tauS,N);
vminA = computeVmin(tauS,dpsiS);

%look finely around vmin from analytics
vin = linspace(vminA-.05,vminA+.05,50);
t  = 0:dt:10;
tV = 0:dt:tmax2;
for i=1:numel(vin)
    [psi,~] = simPsi(JE,vin(i),th,A,C0,psiSamp,tau,N,t,tV,del,0);
    psi = unwrap(psi);
    if psi(end)>dpsiS+dpsiU/2
        moves(i) = 1;
    else
        moves(i) = 0;
    end
end

ii = find(moves<.5,1,'last');
vmin = vin(ii);

end

function linearity = simLinInt(vin,JE,JEmin,JEmax,N,tmax,tmax2,simVars)
A   = simVars.A;
C0  = simVars.C0;
tau = simVars.tau;
dt  = simVars.dt;
del = simVars.del;
th  = simVars.th;
psiSamp = simVars.psiSamp;

[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,dpsiS] = computeDelNact(tauU,tauS,N);

t  = 0:dt:tmax;
tV = 0:dt:tmax2;
for i=1:numel(vin)
    [psi,~] = simPsi(JE,vin(i),th,A,C0,psiSamp,tau,N,t,tV,del,0);
    psi = unwrap(psi);
    [vmax,vmin] = computeBumpVelocity(psi,t,dpsiS,N);
    linearity(i) = vmin./vmax;
end

end

function [evU,evS] = computeEvalsNum(W,Nact)

indsU = 1:Nact;
indsS = 1:Nact-1;

[~,DU] = eig(W(indsU,indsU));
[~,DS] = eig(W(indsS,indsS));

evU = max(real(diag(DU)));
evS = max(real(diag(DS)));

if max(imag(diag(DU)))>1e-10 || max(imag(diag(DS)))>1e-10 
    error('complex evals')
end
end


%----------------------- ANALYTIC TRAJECTORIES ---------------------------%

function [psiS,psiN,vmin] = computePsiStable(t,vin,JE,JEmin,JEmax,tau,N)
% valid for psi<psiN, v<vmin
[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,psiN] = computeDelNact(tauU,tauS,N);
vmin = computeVmin(tauS,psiN);
psiS = vin.*psiN./vmin.*(1-exp(-t./tauS));
end

function [psiU,tU,psiS,tS] = computePsiUnstable(psi0,t,vin,JE,JEmin,JEmax,tau,N)
% valid for psi>psiN, v<vmin
[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,psiN] = computeDelNact(tauU,tauS,N);
vmin = computeVmin(tauS,psiN);

Au = psiN+tauU.*(vmin-vin);
Bu = psi0 - Au;
As = vin.*psiN./vmin;
Bs = psiN.*(1-vin./vmin);

tCross = tauU.*log((psiN-Au)./Bu);
t0 = linspace(0,tCross,1000);
dt = t0(2)-t0(1);
psiU = Au + Bu.*exp(t0./tauU);
t1 = linspace(0,max(t)-t0(end)-dt,1000);
psiS = As + Bs.*exp(-t1./tauS);

tU = t0;
tS = t1+t0(end)+dt;

end

function [psiS,psiU,tS,tU, psiStail,psiUtail,tStail,tUtail] = computePsi(vin,JE,JEmin,JEmax,tau,N)
% find psiS
[~,~,tauU,tauS] = computeTimeConstants(JE,JEmin,JEmax,tau);
[~,psiN] = computeDelNact(tauU,tauS,N);
vmin = computeVmin(tauS,psiN);

tMinS = -tauS.*log(1+vmin/vin); % time when psiS = -psiN
tMaxS = -tauS.*log(1-vmin/vin); % time when psiS = psiN
tS = linspace(tMinS,tMaxS,1000); 
psiS = computePsiStable(tS,vin,JE,JEmin,JEmax,tau,N);

% find psiU
psiThresh = pi/N-psiN;
Au = pi/N - psiThresh.*vin./vmin;
Bu = psiThresh.*(vin./vmin-1);
tauEff = psiThresh./vmin;

tMinU = 0;
tMaxU = tauEff.*log((vin+vmin)./(vin-vmin)); %time when psiU = psiN+2*psiThresh, defined relative to t = 0
tU = linspace(tMinU,tMaxU,1000);
psiU = Au+Bu.*exp(tU./tauEff);

%include tails (for plotting)
f = 20;
tStail   = linspace(tMinS,tMaxS+f*(tMaxS-tMinS),f*1000);
psiStail = computePsiStable(tStail,vin,JE,JEmin,JEmax,tau,N);

tUtail   = linspace(tMinU-f*(tMaxS-tMinS),tMaxU,f*1000);
psiUtail = Au+Bu.*exp(tUtail./tauEff);
end


%---------------------- SIMULATED TRAJECTORIES ---------------------------%

function [psi,h] = simPsi(JE,vdesired,th,A,C0,psiSamp,tau,N,t,tV,del,psi0,addNoise,queNoise)

% if no additive noise was passed, set it to zero
if ~exist('addNoise','var')
    addNoise = 0;
end

% get weight matrices
[Ws,Wa,JI] = getConnectivity(N,JE,A,C0,psiSamp,del);

% if quenched noise was passed, add it to the symmetric weight matrix
if exist('queNoise','var')
    Ws = Ws + queNoise;
end

% initialize bump
h0 = initializeBump(JE,JI,N,th,C0,psi0);

% determine velocity scaling
if vdesired>0
    fslope = getVscale(JE,tV,th,A,C0,psiSamp,del,tau,N); 
    v = vdesired./fslope;
    Wccw = zeros(size(Wa));
    Wcw  = Wa;
else
    v = vdesired;
    Wccw = zeros(size(Wa));
    Wcw  = Wccw;
end

% simulate
[h,~,~,psi,~] = simM(Ws,Wccw,Wcw,C0,tau,-v,h0,t,N,addNoise);

end



function h0 = initializeBump(JE,JI,N,th,C0,psi0)
thcinit = pi/2;
[thc0,~,flag,~] = fzero(@(thc)([0 JE 0]*calculateFunctions([psi0 thc],N)' - 1),thcinit);
if flag < 0
    disp(['Something went wrong, flag for root finder: ' num2str(flag)])
end
[f,~] = calculateFunctions([psi0 thc0],N);
f0    = f(1);
rho0  = -C0/(2*(cos(thc0) + JI*f0));
h0    = 2*rho0*(cos(th-psi0) - cos(thc0));
end

function [vmax,vmin] = computeBumpVelocity(psi,t,dpsiS,N)
dt     = t(2)-t(1);
dpsiU  = pi/N-dpsiS;
[~,i1] = min(abs(psi-dpsiS));
[~,i2] = min(abs(psi-(pi/N+dpsiU)));

dpsi = diff(psi);
vmin = mean([dpsi(i1),dpsi(i1+1)])./dt;
vmax = mean([dpsi(i2),dpsi(i2+1)])./dt;

end

function f = getVscale(JE,t,th,A,C0,psiSamp,del,tau,N,v,addNoise)
%determines scaling for velocity input to match bump velocity

if nargin<11
    addNoise = 0;
    if nargin<10
        v = 50;
    end
end

[Ws,Wa,JI] = getConnectivity(N,JE,A,C0,psiSamp,del);
h0 = initializeBump(JE,JI,N,th,C0,0);

%simulate
[~,~,~,psi,~] = simM(Ws,zeros(size(Wa)),Wa,C0,tau,-v,h0,t,N,addNoise);
psi = unwrap(psi);
%psi = unwrap(psi-psi0);


m = linspace(0,1000,10000);
[~,ii] = min(sum(abs(m.*t' - psi)));
vavg = m(ii);

f = vavg./v;
end

function psiU = unwrap(psi)
% unwraps bump orientation 

thresh = -0.0001;   %threshold for determinining wrapped time points
psiU = psi;         %vector of unwrapped bump orientations
dp   = diff(psiU);
ii   = find(dp<thresh);
for i=1:numel(ii)
    psiU( (ii(i)+1):end ) = psiU( (ii(i)+1):end ) + 2*pi;
end

end
