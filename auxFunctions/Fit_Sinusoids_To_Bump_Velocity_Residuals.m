% This function fits sinusoids of different frequencies to the residual
% bump velocities as a function of bump position

function [BinnedFly] = Fit_Sinusoids_To_Bump_Velocity_Residuals(BinnedFly, Freqs, Bump_HD_Good_Left, Bump_AV_Good_Residual_Left, ...
                            Bump_HD_Good_Right, Bump_AV_Good_Residual_Right, FlyInd)


% Loop over frequencies to fit
A_All_Left=[];   Theta_All_Left=[];  DC_All_Left=[];  Model_FitToData_Left=[];  Model_FitToData_Left_Plot=[];  R2_All_Left=[]; RP_All_Left=[];                    
A_All_Right=[];  Theta_All_Right=[]; DC_All_Right=[]; Model_FitToData_Right=[]; Model_FitToData_Right_Plot=[]; R2_All_Right=[]; RP_All_Right=[];
for fff=1:length(Freqs)
    
    % Update the frequency to fit
    W = Freqs(fff);

    % Define model
    mdl = fittype(@(A, Theta, DC, x) A*sin(W*x + Theta) + DC, 'independent', 'x');
    
    % Fit model to left turns and calculate R^2
    % Contrain amplitudes to be positive to be consistent (i.e. unconstrained, the model could fit a sinsuoid 180 degrees apart but invert the amplitude)
    [fitmdl_left xl] = fit(Bump_HD_Good_Left'/360*2*pi, Bump_AV_Good_Residual_Left', mdl, 'Start', [30, 0, 0],'Lower',[0 -pi -inf], 'Upper',[inf pi inf]);
    Model_FitToData_Left(fff,:) = fitmdl_left(Bump_HD_Good_Left'/360*2*pi);
    Model_FitToData_Left_Plot(fff,:) = fitmdl_left(-pi:pi/1000:pi);
    [R2_Temp_Left P_Temp_Left] = corrcoef(Model_FitToData_Left(fff,:), Bump_AV_Good_Residual_Left);
    R2_Temp_Left=R2_Temp_Left.^2;

    % Fit model to right turns and calculate R^2
    fitmdl_right = fit(Bump_HD_Good_Right'/360*2*pi, Bump_AV_Good_Residual_Right', mdl,'start',[30, 0, 0],'Lower',[0 -pi -inf], 'Upper',[inf pi inf]);
    Model_FitToData_Right(fff,:) = fitmdl_right(Bump_HD_Good_Right'/360*2*pi);
    Model_FitToData_Right_Plot(fff,:) = fitmdl_right(-pi:pi/1000:pi);
    [R2_Temp_Right P_Temp_Right] = corrcoef(Model_FitToData_Right(fff,:), Bump_AV_Good_Residual_Right);
    R2_Temp_Right=R2_Temp_Right.^2;

    % Store parameters
    A_All_Left(fff)       =  fitmdl_left.A;
    Theta_All_Left(fff)   =  fitmdl_left.Theta;
    DC_All_Left(fff)      =  fitmdl_left.DC;
    R2_All_Left(fff)      =  R2_Temp_Left(1,2);
    RP_All_Left(fff)      =  P_Temp_Left(1,2);

    A_All_Right(fff)      =  fitmdl_right.A;
    Theta_All_Right(fff)  =  fitmdl_right.Theta;
    DC_All_Right(fff)     =  fitmdl_right.DC;
    R2_All_Right(fff)     =  R2_Temp_Right(1,2);
    RP_All_Right(fff)     =  P_Temp_Right(1,2);

end


% Store data for all trials
BinnedFly.A_Left(FlyInd,:)        =  A_All_Left;
BinnedFly.A_Right(FlyInd,:)       =  A_All_Right; 
BinnedFly.Theta_Left(FlyInd,:)    =  Theta_All_Left; 
BinnedFly.Theta_Right(FlyInd,:)   =  Theta_All_Right; 
BinnedFly.DC_Left(FlyInd,:)       =  DC_All_Left; 
BinnedFly.DC_Right(FlyInd,:)      =  DC_All_Right;
BinnedFly.R2_Left(FlyInd,:)       =  R2_All_Left;
BinnedFly.R2_Right(FlyInd,:)      =  R2_All_Right;
BinnedFly.RP_Left(FlyInd,:)       =  RP_All_Left;
BinnedFly.RP_Right(FlyInd,:)      =  RP_All_Right;
BinnedFly.Model_FitToData_Left_Plot{FlyInd}  = Model_FitToData_Left_Plot;
BinnedFly.Model_FitToData_Right_Plot{FlyInd} = Model_FitToData_Right_Plot;








