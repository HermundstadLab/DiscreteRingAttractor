% This function fits sinusoids of different frequencies to the drift at
% different bump positions (starting and restarting)

function [A_All_Begin, A_All_End, Theta_All_Begin, Theta_All_End, DC_All_Begin, DC_All_End,...
          Model_FitToData_Begin, Model_FitToData_End, Model_FitToData_Begin_Plot, Model_FitToData_End_Plot,...
          R2_All_Begin, R2_All_End, RP_All_Begin, RP_All_End] = Fit_Sinusoids_DriftHists(Freqs, ...
          EBPositions, Bump_StopBegin, Bump_StopEnd)


% Fit sinusoids of different frequencies to the HD x residual bump velocity plots, conditioned on turn direction.
A_All_Begin=[];                 A_All_End=[];
Theta_All_Begin=[];             Theta_All_End=[];
DC_All_Begin=[];                DC_All_End=[];
Model_FitToData_Begin=[];       Model_FitToData_End=[];
Model_FitToData_Begin_Plot=[];  Model_FitToData_End_Plot=[];
R2_All_Begin=[];                R2_All_End=[];
RP_All_Begin=[];                RP_All_End=[];
for fff=1:length(Freqs)
    
    % Update the frequency to fit
    W = Freqs(fff);

    % Define model
    mdl = fittype(@(A, Theta, DC, x) A*sin(W*x + Theta) + DC, 'independent', 'x');

    % % Shift BumpStopBegin and BumpStopEnd, both probability density
    % % functions, by their mean, for a regular sine fit (rather than constraining
    % % sine to be positive at all positions.
    inclStart = ~isnan(Bump_StopBegin);
    exclStart = isnan(Bump_StopBegin);
    inclStop = ~isnan(Bump_StopEnd);
    exclStop = isnan(Bump_StopEnd);
    Bump_StopBegin = Bump_StopBegin(inclStart);
    Bump_StopBegin = Bump_StopBegin-mean(Bump_StopBegin);
    Bump_StopEnd = Bump_StopEnd(inclStop);
    Bump_StopEnd = Bump_StopEnd-mean(Bump_StopEnd);

    % Fit model to positions at the beginning of stops and calculate R^2
    % Constrain amplitudes to be positive to be consistent (i.e. unconstrained, the model could fit a sinsuoid 180 degrees apart but invert the amplitude)

    EBPositionsStart = EBPositions';
    EBPositionsStart = EBPositionsStart(inclStart);

    [fitmdl_begin, gof_begin] = fit(EBPositionsStart, Bump_StopBegin', mdl, 'Normalize', 'on', 'MaxIter', 10000, 'Start', [0.0001, 0, 0],'Lower',[0 -pi -2], 'Upper',[1 pi 2]);
    Model_FitToData_Begin(fff,:) = fitmdl_begin(EBPositionsStart);
    Model_FitToData_Begin_Plot(fff,:) = fitmdl_begin(-pi:pi/1000:pi);
    R2_Temp_Begin=gof_begin.rsquare;
    [~, P_Temp_Begin] = corrcoef(Model_FitToData_Begin(fff,:), Bump_StopBegin);

    % Fit model to End turns and calculate R^2
    EBPositionsStop = EBPositions';
    EBPositionsStop = EBPositionsStop(inclStop);
    [fitmdl_end, gof_end] = fit(EBPositionsStop, Bump_StopEnd', mdl, 'Normalize', 'on', 'MaxIter', 10000, 'Start',[0.0001, 0, 0],'Lower',[0 -pi -2], 'Upper',[1 pi 2]);
    Model_FitToData_End(fff,:) = fitmdl_end(EBPositionsStop);
    Model_FitToData_End_Plot(fff,:) = fitmdl_end(-pi:pi/1000:pi);
    [~, P_Temp_End] = corrcoef(Model_FitToData_End(fff,:), Bump_StopEnd);
    R2_Temp_End=gof_end.rsquare;

    % Store parameters
    A_All_Begin(fff)       =  fitmdl_begin.A;
    Theta_All_Begin(fff)   =  fitmdl_begin.Theta;
    DC_All_Begin(fff)      =  fitmdl_begin.DC;
    R2_All_Begin(fff)      =  R2_Temp_Begin;
    RP_All_Begin(fff)      =  P_Temp_Begin(1,2);

    A_All_End(fff)      =  fitmdl_end.A;
    Theta_All_End(fff)  =  fitmdl_end.Theta;
    DC_All_End(fff)     =  fitmdl_end.DC;
    R2_All_End(fff)     =  R2_Temp_End;
    RP_All_End(fff)     =  P_Temp_End(1,2);

end