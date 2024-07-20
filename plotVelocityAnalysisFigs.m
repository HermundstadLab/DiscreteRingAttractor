function plotVelocityAnalysisFigs(DarkData)

%% This script reproduces the bump velocity results in Figure 1 of Noorman et al 2024

% Concatonate each fly's data across trials
[FlyData]=Combine_Data_Across_Trials(DarkData);

% Make directory to save plots to
PlotDir = 'velocityAnalysis_Figures/';
if ~isdir(PlotDir); mkdir(PlotDir); end


%% Step 1: Compute the residual bump velocity for each fly's data.

% Configurable variables
vlow  = 5;       % Lower threshold on angular velocity (used to subset data; vlow = 5 for paper)
vhigh = 300;     % Upper threshold on angular velocity (used to subset data; vhigh = 300 for paper) 

% Get sample period (data was sampled at 10Hz, so dT=0.1 s).
dT = 1/FlyData.TwoP_Fs(1);

% Loop over flies and analyze the data
AllFly.Bump_HD_Good_Left={};            AllFly.Bump_HD_Good_Right={};
AllFly.Bump_AV_Good_Residual_Left={};   AllFly.Bump_AV_Good_Residual_Right={};
for TrialInd=1:size(FlyData.BumpR,1)
    
    display(['Fly ' num2str(TrialInd) ' of ' num2str(size(FlyData.BumpR,1))])

    %%%%%%%%%%%%%%%% Get fly/bump velocities during turns %%%%%%%%%%%%%%%%%
    % For this analysis, the bump velocity (Bump_AV) can only be computed
    % during periods where there is a clear activity bump, so we'll use a
    % subset of the data when the bump's mean resultant vector length is
    % large (BumpR>0.125) and when the fly is turning (from 5 to 300 deg/s)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get fly's rotational velocity and compute the fly's HD from it
    FlyAV  = FlyData.FlyRotate_w2p(TrialInd,:);    % in deg/s     

    % Unwrap the bump's position and compute the bump velocity
    BumpHD = FlyData.BumpPosition(TrialInd,:)/2/pi*360;       % in deg 
    BumpHD_Unwrap =  unwrap(BumpHD/360*2*pi, pi)/2/pi*360;    % in deg
    Bump_AV = -diff(BumpHD_Unwrap)/dT;                        % in deg/s
    Bump_AV = mean([Bump_AV 0; 0 Bump_AV], 1); 
    
    % Subset the fly/bump AV estimates by bump mean resultant vector length and the fly's AV
    BumpR = FlyData.BumpR(TrialInd,:);
    Logica = (BumpR>0.125) & (abs(FlyAV) > vlow & abs(FlyAV) < vhigh);
    FlyAV_Good = FlyAV(Logica);
    Bump_AV_Good = Bump_AV(Logica);
    Bump_HD_Good = BumpHD(Logica);


    %%%%%%%%%%%%%% Compute naive gain for left/right turns %%%%%%%%%%%%%%%%
    % Each fly has a its own naive gain between rotations on the ball and
    % rotations of its bump. To account for this, we need to perform a
    % linear regression between fly velocity and bump velocity, and the
    % residual of the linear fit indicates whether the bump is moving
    % faster or slower than expected. We do this for left turns and right
    % turns separately. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get the fly's AV and Bump's AV/HD for left turns
    LeftTurnLogica      =  FlyAV_Good<0;
    FlyAV_Good_Left     =  FlyAV_Good(LeftTurnLogica);
    Bump_AV_Good_Left   =  Bump_AV_Good(LeftTurnLogica);
    Bump_HD_Good_Left   =  Bump_HD_Good(LeftTurnLogica);

    % Perform a linear regression (OLS) between the fly's AV and bump's AV for left turns
    LinearFit_Left      =  polyfit(FlyAV_Good_Left, Bump_AV_Good_Left, 1);

    % Compute the residuals
    Bump_AV_Good_Residual_Left  =  Bump_AV_Good_Left - (FlyAV_Good_Left*LinearFit_Left(1) + LinearFit_Left(2));
   
    % Invert so positive residual is bump going faster than fly (for left turns)
    Bump_AV_Good_Residual_Left  = -Bump_AV_Good_Residual_Left; 

    % Get the fly's AV and Bump's AV/HD for right turns
    RightTurnLogica     =  FlyAV_Good>0; 
    FlyAV_Good_Right    =  FlyAV_Good(RightTurnLogica);
    Bump_AV_Good_Right  =  Bump_AV_Good(RightTurnLogica);
    Bump_HD_Good_Right  =  Bump_HD_Good(RightTurnLogica);

    % Perform a linear regression (OLS) between the fly's AV and bump's AV for right turns
    LinearFit_Right     =  polyfit(FlyAV_Good_Right, Bump_AV_Good_Right, 1);

    % Compute the residuals
    Bump_AV_Good_Residual_Right =  Bump_AV_Good_Right - (FlyAV_Good_Right*LinearFit_Right(1) + LinearFit_Right(2));

    % Plot the linear fits for left/right turns
    Plot_Fly_Vs_Bump_Velocity(FlyAV_Good, Bump_AV_Good, LinearFit_Right, LinearFit_Left, TrialInd, PlotDir);

    % Add this fly's data to the AllFly data structure
    AllFly.Bump_HD_Good_Left{TrialInd}            =  Bump_HD_Good_Left;
    AllFly.Bump_HD_Good_Right{TrialInd}           =  Bump_HD_Good_Right;
    AllFly.Bump_AV_Good_Residual_Left{TrialInd}   =  Bump_AV_Good_Residual_Left;
    AllFly.Bump_AV_Good_Residual_Right{TrialInd}  =  Bump_AV_Good_Residual_Right;

end



%% Step 2: Bin residual bump velocity by bump position to get an estimate of
% whether the bump moves faster/slower than expected, on average, at each
% bump position.

% Bin bump position and residuals into bins of 5.625 degrees (i.e. 360/64)
Bins           =  64;    % number of HD bins (64 for paper)
HD_Bin_Edges   =  -180:(360/Bins):180;
HD_Bin_Centers =  HD_Bin_Edges(1:end-1) + 360/Bins/2;
[Binned] = Bin_Bump_Velocity_Residuals(AllFly, HD_Bin_Edges);



%% Step 3: Fit 8 and 16 Hz sinusoids to the binned residual bump velocity as 
%  a function of bump position.

% Frequencies of sinusoids to fit
Freqs=[8, 16];

% Loop over flies and fit sinusoids of 8 and 18 Hz
BinnedFly.A_Left=[];  BinnedFly.Theta_Left=[];  BinnedFly.DC_Left=[];  BinnedFly.R2_Left=[];  BinnedFly.RP_Left=[];   
BinnedFly.A_Right=[]; BinnedFly.Theta_Right=[]; BinnedFly.DC_Right=[]; BinnedFly.R2_Right=[]; BinnedFly.RP_Right=[];  
for FlyInd = 1 : length(AllFly.Bump_HD_Good_Left)

    % Get this fly's binned bump position and residual bump velocities 
    Bump_HD_Good_Left            =  Binned.Bump_HD_Good_Left(FlyInd,:);
    Bump_AV_Good_Residual_Left   =  Binned.Bump_AV_Good_Residual_Left(FlyInd,:);
    Bump_HD_Good_Right           =  Binned.Bump_HD_Good_Right(FlyInd,:);
    Bump_AV_Good_Residual_Right  =  Binned.Bump_AV_Good_Residual_Right(FlyInd,:);

    % Fit sinusoids to the residual bump velocity as a function of bump position
    [BinnedFly] = Fit_Sinusoids_To_Bump_Velocity_Residuals(BinnedFly, Freqs, Bump_HD_Good_Left, ...
        Bump_AV_Good_Residual_Left, Bump_HD_Good_Right, Bump_AV_Good_Residual_Right, FlyInd);
                        
end



%% Step 4: Plot the results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Plot the results for Figure 1J %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

subplot(2,1,1)
hold on;

plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Left, 'color', [0.5 0.5 0.5])
h1 = plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Left(2,:), 'color', [1 0 0]);
h2 = plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Left(6,:), 'color', [0 0 1]);
plot(HD_Bin_Centers, mean(Binned.Bump_AV_Good_Residual_Left,1), 'k', 'linewidth', 2)
ylim([-30 30])
xlim([-180 180])
xlabel('bump position (deg)')
ylabel('residual bump velocity (deg/s)')
title('Left Turns')
set(gca, 'xtick', [-180 -90 0 90 180], 'ytick', [-30 -15 0 15 30])
legend([h1 h2],'FlyA','FlyB');

subplot(2,1,2)
hold on;
plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Right, 'color', [0.5 0.5 0.5])
h1 = plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Right(2,:), 'color', [1 0 0]);
h2 = plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Right(6,:), 'color', [0 0 1]);
plot(HD_Bin_Centers, mean(Binned.Bump_AV_Good_Residual_Right,1), 'k', 'linewidth', 2)
ylim([-30 30])
xlim([-180 180])
xlabel('bump position (deg)')
ylabel('residual bump velocity (deg/s)')
title('Right Turns')
set(gca, 'xtick', [-180 -90 0 90 180], 'ytick', [-30 -15 0 15 30])
legend([h1 h2],'FlyA','FlyB');

set(gcf,'PaperOrientation','portrait','PaperPosition',[0 0 5 10],'PaperSize',[5 10]) % Vectors are width then height
print(gcf, '-dpdf', '-r350', [PlotDir '\Figure_1J_BinnedResiduals'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% R^2 value of fits reported Figure S2 B  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 8 Hz stats
EightHz.Ind = find(Freqs == 8);
EightHz.R2_Left  = BinnedFly.R2_Left(:, EightHz.Ind);
EightHz.RP_Left  = BinnedFly.RP_Left(:, EightHz.Ind);
EightHz.R2_Right = BinnedFly.R2_Right(:, EightHz.Ind);
EightHz.RP_Right = BinnedFly.RP_Right(:, EightHz.Ind);


% 16 Hz stats
SixteenHz.Ind = find(Freqs == 16);
SixteenHz.R2_Left  = BinnedFly.R2_Left(:, SixteenHz.Ind);
SixteenHz.RP_Left  = BinnedFly.RP_Left(:, SixteenHz.Ind);
SixteenHz.R2_Right = BinnedFly.R2_Right(:, SixteenHz.Ind);
SixteenHz.RP_Right = BinnedFly.RP_Right(:, SixteenHz.Ind);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Plot results in Figure S2  A %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure;
ModelFitBins= (-pi:pi/1000:pi)/2/pi*360; 
for fff=1:10;

    % Plot left turn data
    subplot(10,2,fff*2 -1);
    hold on;
    plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Left(fff,:), 'color', [0.5 0.5 0.5])
    plot(ModelFitBins,   BinnedFly.Model_FitToData_Left_Plot{fff}(EightHz.Ind,:), 'color', 'r')
    plot(ModelFitBins,   BinnedFly.Model_FitToData_Left_Plot{fff}(SixteenHz.Ind,:), 'color', 'b')
    xlim([-180 180])
    ylim([-30 30])
    set(gca, 'xtick', [-180:90:180], 'ytick', [-30:15:30])
    title(['Fly ' num2str(fff) ' Left Turns'])
    xlabel('bump orientation (deg)')
    ylabel({'residual bump'; 'velocity (deg/s)'})

    % Plot right turn data
    subplot(10,2,fff*2);
    hold on;
    plot(HD_Bin_Centers, Binned.Bump_AV_Good_Residual_Right(fff,:), 'color', [0.5 0.5 0.5])
    plot(ModelFitBins,   BinnedFly.Model_FitToData_Right_Plot{fff}(EightHz.Ind,:), 'color', 'r')
    plot(ModelFitBins,   BinnedFly.Model_FitToData_Right_Plot{fff}(SixteenHz.Ind,:), 'color', 'b')
    xlim([-180 180])
    ylim([-30 30])
    set(gca, 'xtick', [-180:90:180], 'ytick', [-30:15:30])
    title(['Fly ' num2str(fff) ' Right Turns'])
    xlabel('bump orientation (deg)')
    ylabel({'residual bump'; 'velocity (deg/s)'})

end

set(gcf,'PaperOrientation','portrait','PaperPosition',[0 0 8 16],'PaperSize',[8 16]) % Vectors are width then height
print(gcf, '-dpdf', '-r350', [PlotDir '\Figure_S2_A_BumpResiduals'])





