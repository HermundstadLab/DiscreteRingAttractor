%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN SCRIPT                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------ KEY --------------------------------------%
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

%to generate all figures, run: plotMSfigs(1:27)
%-------------------------------------------------------------------------%


%% FIGURE 2A 
plotMSfigs(1);        %plot connectivity matrix schematic

%% FIGURE 2B
plotMSfigs(2);        %plot optimal JE          

%% FIGURE 2E
plotMSfigs(3);        %plot energy heatmaps 

%% Figure 2F-G
plotMSfigs([4,5]);    %plot bump trajectories with and without input   

%% FIGURE 2H
plotMSfigs(6);        %plot connectivity matrices for different JE              

%% FIGURE 3C
plotMSfigs(7);        %plot bump trajectories in different regimes, without input           

%% FIGURE 3D
plotMSfigs(8);        %plot drift rates in different regimes                

%% FIGURE 3E-F
plotMSfigs(9);        %plot simplified energy landscape, without input               

%% FIGURE 3G, ED FIGURE 8G
plotMSfigs(10);       %plot bump trajectories in different regimes, varying JE, without input               

%% FIGURE 3H
plotMSfigs(11);       %plot net drift speed             

%% FIGURE I-J
plotMSfigs(12);       %plot simplified energy landscape, with input              

%% FIGURE 3K, ED FIGURE 8G
plotMSfigs(13);       %plot bump trajectories in different regimes, varying JE, with input  

%% FIGURE 3L
plotMSfigs(14);  %plot threshold velocity                        

%% FIGURE 3L
plotMSfigs(15);       %plot linearity of integration                 

%% FIGURE 4A-C, ED FIGURE 9A
plotMSfigs(16);       %plot robustness to parameter tuning                  

%% FIGURE 4D-F, ED FIGURE 9B
plotMSfigs(17);       %plot noise robustness                

%% ED FIGURE 3A
plotMSfigs(18);       %plot phase diagram                 

%% ED FIGURE 3B-C
plotMSfigs(19);       %plot fixed point conditions                

%% ED FIGURE 3D-I
plotMSfigs(20);       %plot linearization  

%% ED FIGURE 4
plotMSfigs(21);       %plot eigenvalues of Hessian

%% ED FIGURE 5
run plotEDfig5.m

%% ED FIGURE 6
plotMSfigs(22);       %plot leading eigenvalues without input

%% ED FIGURE 7
plotMSfigs(23);       %plot leading eigenvalues with small input

%% ED FIGURE 8B
plotMSfigs(24);       %plot drift speeds in different regimes, grayscale

%% ED FIGURE 8C,E
plotMSfigs(25);       %plot example fixed points in different regimes as a function of input velocity 

%% ED FIGURE 8D
plotMSfigs(26);       %plot all fixed points in different regimes for velocities below threshold

%% ED FIGURE 8F
plotMSfigs(27);       %plot all fixed points in different regimes for velocities above threshold
