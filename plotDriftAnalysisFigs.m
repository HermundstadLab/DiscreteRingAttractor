function [analysisResults] = plotDataAnalysisFigs(DarkData)
% Perform drift analyses associated with Noorman et al. 2024 figures.

% First, make a couple of sample plots of bump movement.
% Examples used in Figure 1 of the paper below:
% Fly 378, Trial 5, Segment 26
Plot_Bump_Examples(DarkData,10,26);
% Fly 390, Trial 4, Segment 4
Plot_Bump_Examples(DarkData,29,4);

allFlyIdxs = unique(DarkData.Fly);
nTotalFlies = length(allFlyIdxs);
nTotTrials = length(DarkData.Fly);

% Set analysis parameters
tFrame = 0.1; % seconds/frame
threshNoRot = eps; 
threshNoFwd = eps;
tStill = 0.3; % seconds with 0 ang vel to be considered not rotating
EBsubdiv = 64;

nRows = 2;
nCols = 5;

% Initialize output struct
for k = 1:nTotalFlies
    analysisResults(k).diffBumpPos = [];
    analysisResults(k).diffTime = [];
    analysisResults(k).bumpPos = [];
    analysisResults(k).restartBumpPos = [];
    analysisResults(k).beginIdxs = [];
    analysisResults(k).endIdxs = [];

    analysisResults(k).RestartModelFit = [];
    analysisResults(k).RestartModelFitR2 = [];

    analysisResults(k).StartDriftModelFit = [];
    analysisResults(k).StartDriftModelFitR2 = [];
    analysisResults(k).StartDriftModelFitPval = [];
    analysisResults(k).RestartDriftModelFit = [];
    analysisResults(k).RestartDriftModelFitR2 = [];
    analysisResults(k).RestartDriftModelFitPval = [];

    analysisResults(k).watsonsU2TestPVal = [];
    analysisResults(k).watsonsU2TestU2Obs = [];
    analysisResults(k).watsonsU2TestU2H0W = [];
end

currFlyNum = 0;
currFlyIdx = 0;
count = 0;

% Go through all trials, picking out the 5 that correspond to each fly
for i = 1:nTotTrials
    potentialFlyIdx = DarkData.Fly(i);

    if potentialFlyIdx ~= currFlyIdx
        currFlyNum = currFlyNum + 1;
        currFlyIdx = potentialFlyIdx;
        nTrialsThisFly = 0;
        count = 0;
    end

    nTrialsThisFly = nTrialsThisFly + 1;

    dfwd = DarkData.FlyForward_w2p(i,:);
    drot = DarkData.FlyRotate_w2p(i,:);

    % Set all below thresh values to 0 or 1 depending on speed.
    % This will allow us to extract the beginnings and ends of 
    % periods of standing still.
    dfwd(abs(dfwd) <= threshNoFwd) = 0;
    dfwd(abs(dfwd) > threshNoFwd) = 1;

    drot(abs(drot) <= threshNoRot) = 0;
    drot(abs(drot) > threshNoRot) = 1;
    
    % Handle edge cases at the beginning of the trial
    if (drot(1) == 0)
        % fly's not rotating at the beginning of the trial, so consider
        % that the potential beginning of the first bout of standing still
        beginIdxs = [1 find(drot(1:end-1) - drot(2:end) > 0)];
    else
        beginIdxs = find(drot(1:end-1) - drot(2:end) > 0);
    end
    
    % Handle edge cases at the end of the trial
    if (drot(end) == 0)
        % Fly is not rotating at the end of the trial, so consider it the
        % potential end of the last bout of standing still
        endIdxs = [find(drot(2:end) - drot(1:end-1) > 0) length(drot)];
    else
        endIdxs = find(drot(2:end) - drot(1:end-1) > 0);
    end

    bumpPosTmp = DarkData.BumpPosition(i,:);

    % Extract beginning and end bump positions for all standing bouts
    analysisResults(currFlyNum).beginIdxs = [];
    analysisResults(currFlyNum).endIdxs = [];
    if (length(beginIdxs) == length(endIdxs))
        beginBumpPos = bumpPosTmp(beginIdxs);
        endBumpPos = bumpPosTmp(endIdxs);
        for j = 1:length(endIdxs)
            % We know the fly isn't turning during this epoch, but let's
            % check to make sure that fly isn't walking forward either.
            if isempty(find(dfwd(beginIdxs(j):endIdxs(j)),1))
                count = count + 1;
                analysisResults(currFlyNum).beginIdxs = [analysisResults(currFlyNum).beginIdxs beginIdxs(j)];
                analysisResults(currFlyNum).endIdxs = [analysisResults(currFlyNum).endIdxs endIdxs(j)];
                analysisResults(currFlyNum).diffBumpPos(count) = circ_dist(endBumpPos(j), beginBumpPos(j));

                analysisResults(currFlyNum).diffTime(count) = tFrame*(endIdxs(j) - beginIdxs(j));
                analysisResults(currFlyNum).bumpPos(count) = DarkData.BumpPosition(i, beginIdxs(j));
                analysisResults(currFlyNum).restartBumpPos(count) = DarkData.BumpPosition(i, endIdxs(j));
            end
        end
    else
        disp(["idx error in Trial:" i]);
    end
end

dTimeAll = [];
dBumpPosAll = [];

figure;
% Plot drift for different durations of standing bouts
for k = 1:nTotalFlies
    subplot(nRows+1,nCols,k);

    dTime = analysisResults(k).diffTime;
    dTimeAll = [dTimeAll dTime];

    dBumpPos = analysisResults(k).diffBumpPos;
    dBumpPosAll = [dBumpPosAll dBumpPos];
    scatter(dTime(dTime>0.2), dBumpPos(dTime>0.2), 'filled');
    ylim([-pi pi]);
    yticks([-pi -pi/2 0 pi/2 pi]);
    xlabel('Time (s)');
    ylabel('drift (rad)');
    title(['Fly #' num2str(k)]);
end

% Plot scatter of drift across flies
subplot(nRows+1,nCols,k+1);
scatter(dTimeAll, dBumpPosAll, 'filled');
yticks([-pi -pi/2 0 pi/2 pi]);
xlabel('Time (s)');
ylabel('drift (rad)');
title('All flies drift vs time standing');
ylim([-pi pi]);

dTimeAll = [];
dBumpPosAll = [];

figure;
% Plot histogram of extent of drifts 
for k = 1:nTotalFlies
    subplot(nRows+1,nCols,k);

    dTime = analysisResults(k).diffTime;
    dBumpPos = analysisResults(k).diffBumpPos;

    histogram(dBumpPos(dTime > tStill), -pi:2*pi/61:pi);
    xlim([-pi pi]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xlabel('Drift (rad)');
    ylabel('number of epochs');
    title(['Fly #' num2str(k) '; Tstanding > ' num2str(tStill) 's']);
    dTimeAll = [dTimeAll dTime];
    dBumpPosAll = [dBumpPosAll dBumpPos];
end

% Plot drift histograms across flies
subplot(nRows+1,nCols,k+1);
histogram(dBumpPosAll(dTimeAll > tStill), -pi:2*pi/61:pi);
xlim([-pi pi]);
xlabel('Drift (rad)');
ylabel('number of epochs');
xticks([-pi -pi/2 0 pi/2 pi]);
title(['All flies; Tstanding > ' num2str(tStill) 's']);

% Now plot magnitude of drift as a function of bump position in the EB.
% Select bouts of stopping ranging between min and max times that ensure
% sufficient sampling of each position and consistency across flies.
figure;
tStillMax = 2; % max seconds with 0 ang vel to be included in computation of EB-position-dependent drift
EBEdges = -pi:pi/EBsubdiv:pi; % Bins for position-dependent drift
EBPositions = EBEdges+pi/(EBsubdiv*2);
EBPositions = EBPositions(1:end-1);
nBoutsStandingAll = 0;
sumDriftStart = zeros(1, length(EBEdges)-1);
binFlyCount = zeros(1, length(EBEdges)-1);
for k = 1:nTotalFlies
    subplot(nRows,nCols,k);
    dTime = analysisResults(k).diffTime;
    dBumpPos = analysisResults(k).diffBumpPos;
    % We'll bin using EB positions where the bump is when the fly begins
    % a standing bout, as well as (separately) the bump's position in the 
    % EB when the fly resumes walking at the end of a standing bout.
    bumpPos = analysisResults(k).bumpPos;
    bumpPosRestart = analysisResults(k).restartBumpPos;
    
    threshIdxs = find(dTime > tStill & dTime < tStillMax);
    dBumpPos = dBumpPos(threshIdxs);
    dTime = dTime(threshIdxs);
    bumpPos = bumpPos(threshIdxs);
    bumpPosRestart = bumpPosRestart(threshIdxs);
    nBoutsStanding = length(threshIdxs);
    nBoutsStandingAll = nBoutsStandingAll + nBoutsStanding;

    scatter(bumpPos, dBumpPos, 'filled');
    hold on;

    [~, ~, binsStart] = histcounts(bumpPos, EBEdges);
    driftStart = nan*ones(1, length(EBEdges)-1);
    for i = 1:length(EBEdges)-1
        dStart_this = dBumpPos(binsStart == i);
        if ~isnan(dStart_this)
            driftStart(i) = mean(dStart_this);
            sumDriftStart(i) = sumDriftStart(i) + driftStart(i);
            binFlyCount(i) = binFlyCount(i)+1;
        end
    end

    [~, ~, binsRestart] = histcounts(bumpPosRestart, EBEdges);
    driftRestart = nan*ones(1, length(EBEdges)-1);
    for i = 1:length(EBEdges)-1
        dRestart_this = dBumpPos(binsRestart == i);
        if ~isnan(dRestart_this)
            driftRestart(i) = mean(dRestart_this);
        end
    end

    Freqs = [8 16];
    [A_All_Begin, A_All_End, Theta_All_Begin, Theta_All_End, DC_All_Begin, DC_All_End,...
        Model_FitToData_Begin, Model_FitToData_End, Model_FitToData_Begin_Plot, Model_FitToData_End_Plot,...
        R2_All_Begin, R2_All_End, RP_All_Begin, RP_All_End] = ...
        Fit_Sinusoids_DriftHists(Freqs, EBPositions(binsStart), bumpPos(binsStart),...
        bumpPos(binsStart)); % dummy input only to match function

    analysisResults(k).StartDriftModelFit = Model_FitToData_Begin;
    analysisResults(k).StartDriftModelFitR2 = R2_All_Begin;
    analysisResults(k).StartDriftModelFitPval = RP_All_Begin;
    analysisResults(k).RestartDriftModelFit = Model_FitToData_End;
    analysisResults(k).RestartDriftModelFitR2 = R2_All_End;
    analysisResults(k).RestartDriftModelFitPval = RP_All_End;

    plot(EBPositions, driftStart, 'b', 'LineWidth', 2);
    hold on;
    plot(-pi:pi/1000:pi, Model_FitToData_Begin_Plot');
    ylim([-pi pi]);
    xlim([-pi pi]);
    xlabel('EB pos (rad)');
    ylabel('drift (rad)');
    txt8 = ['8 R2 = ' num2str(R2_All_Begin(1))];
    txt16 = ['16 R2 = ' num2str(R2_All_Begin(2))];
    legend({'data','mean',txt8, txt16});
    title(['Fly #' num2str(k) '; Tstanding > ' num2str(tStill) 's']);
end

fBegin = figure;
fEnd = figure;
fCDFs = figure;
for k = 1:nTotalFlies
    beginBumpPosSet = analysisResults(k).bumpPos;
    endBumpPosSet = analysisResults(k).restartBumpPos;

    figure(fCDFs);
    subplot(nRows,nCols,k);
    cdfplot(beginBumpPosSet);
    hold on;
    cdfplot(endBumpPosSet);
    legend('begin pos', 'end pos');
    xlabel('Bump position (rads)');
    ylabel('F(x)');

    NPerms = 500;
    [pvalW2,U2_obsW,U2_H0W]=watsons_U2_perm_test(beginBumpPosSet,endBumpPosSet,NPerms);

    analysisResults(k).watsonsU2TestPVal = pvalW2;
    analysisResults(k).watsonsU2TestU2Obs = U2_obsW;
    analysisResults(k).watsonsU2TestU2H0W = U2_H0W;

    figure(fBegin);
    subplot(nRows,nCols,k);
    histogram(beginBumpPosSet, EBsubdiv);
    xlabel('Begin stop bump position');
    ylabel('Binned frequency');

    figure(fEnd);
    subplot(nRows,nCols,k);
    histogram(endBumpPosSet, EBsubdiv);
    xlabel('End stop bump position');
    ylabel('Binned frequency');
end

% Now plot position of bump when it restarts as a function of EB position
bumpPosFig = figure;
EBEdges = -pi:pi/EBsubdiv:pi;
sumRestart = zeros(1, length(EBEdges)-1);
binFlyCount = zeros(1, length(EBEdges)-1);
nBoutsStandingAll = 0;

for k = 1:nTotalFlies
    subplot(nRows,nCols,k);

    dTime = analysisResults(k).diffTime;
    bumpPos = analysisResults(k).bumpPos;
    restartBumpPos = analysisResults(k).restartBumpPos;
    threshIdxs = find(dTime > tStill);
    bumpPos = bumpPos(threshIdxs);
    restartBumpPos = restartBumpPos(threshIdxs);
    nBoutsStanding = length(threshIdxs);
    nBoutsStandingAll = nBoutsStandingAll + nBoutsStanding;

    figure(bumpPosFig);
    scatter(bumpPos, restartBumpPos, 'filled');
    hold on;

    [~, ~, bins] = histcounts(bumpPos, EBEdges);
    restartPos = nan*ones(1, length(EBEdges)-1);
    for i = 1:length(EBEdges)-1
        restart_this = restartBumpPos(bins==i);
        if ~isnan(restart_this)
            restartPos(i) = circ_mean(restart_this');
            sumRestart(i) = sumRestart(i) + restartPos(i);
            binFlyCount(i) = binFlyCount(i)+1;
        end
    end

    axis equal;
    ylim([-pi pi]);
    xlim([-pi pi]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    yticks([-pi -pi/2 0 pi/2 pi]);
    xlabel('EB pos (rad)');
    ylabel('EB pos at restart (rad)');
    title(['Fly#' num2str(k) '; Tstand > ' num2str(tStill) 's; nBouts = ' num2str(nBoutsStanding)]);
end

