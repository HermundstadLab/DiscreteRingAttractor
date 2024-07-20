%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN SCRIPT                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To generate all figure panels from the paper, run the following cells in order. 
% Please refer to the key at bottom of this script for plot indices corresponding
% to individual figure panels called by plotModelingFigs(index)

%% SET THE DIRECTORY WHERE DATA IS STORED
dataDir = 'data/';

%% ADD DIRECTORIES AND SUBDIRECTORIES TO CURRENT PATH
addpath("alternateNetworkSimulations/","auxFunctions/","noiseSimulations/",...
  "auxFunctions/perceptually_uniform_colormaps_v1.3.2/Colormaps (5)/Colormaps",...
  "auxFunctions/redblue_v1.0.1","auxFunctions/brewermap","auxFunctions/watsons_U2")

%% FIGURE 1E,H-J; ED FIGS 1-2
data = load([dataDir,'NoormanEtAlFlyWalkingInDarkness2PBehData.mat']);
DarkData = data.DarkData;
[analysisResults] = plotDriftAnalysisFigs(DarkData);

%% FIGURE 2A 
plotModelingFigs(1);        %plot connectivity matrix schematic

%% FIGURE 2D
plotModelingFigs(2);        %plot optimal JE          

%% FIGURE 2E
plotModelingFigs(3);        %plot energy heatmaps 

%% Figure 2F-G
plotModelingFigs([4,5]);    %plot bump trajectories with and without input   

%% FIGURE 2H
plotModelingFigs(6);        %plot connectivity matrices for different JE              

%% FIGURE 3C
plotModelingFigs(7);        %plot bump trajectories in different regimes, without input           

%% FIGURE 3D
plotModelingFigs(8);        %plot drift rates in different regimes                

%% FIGURE 3E-F
plotModelingFigs(9);        %plot simplified energy landscape, without input               

%% FIGURE 3G, ED FIGURE 8G
plotModelingFigs(10);       %plot bump trajectories in different regimes, varying JE, without input               

%% FIGURE 3H
plotModelingFigs(11);       %plot net drift speed             

%% FIGURE I-J
plotModelingFigs(12);       %plot simplified energy landscape, with input              

%% FIGURE 3K, ED FIGURE 8G
plotModelingFigs(13);       %plot bump trajectories in different regimes, varying JE, with input  

%% FIGURE 3L
plotModelingFigs(14);  %plot threshold velocity                        

%% FIGURE 3L
plotModelingFigs(15);       %plot linearity of integration                 

%% FIGURE 4A-C, ED FIGURE 9A
plotModelingFigs(16);       %plot robustness to parameter tuning                  

%% FIGURE 4D-F, ED FIGURE 9B
plotModelingFigs(17);       %plot noise robustness                

%% ED FIGURE 3A
plotModelingFigs(18);       %plot phase diagram                 

%% ED FIGURE 3B-C
plotModelingFigs(19);       %plot fixed point conditions                

%% ED FIGURE 3D-I
plotModelingFigs(20);       %plot linearization  

%% ED FIGURE 4
plotModelingFigs(21);       %plot eigenvalues of Hessian

%% ED FIGURE 5
run plotAltNetworkSims.m    %run simulations for alternate networks and plot results

%% ED FIGURE 6
plotModelingFigs(22);       %plot leading eigenvalues without input

%% ED FIGURE 7
plotModelingFigs(23);       %plot leading eigenvalues with small input

%% ED FIGURE 8B
plotModelingFigs(24);       %plot drift speeds in different regimes, grayscale

%% ED FIGURE 8C,E
plotModelingFigs(25);       %plot example fixed points in different regimes as a function of input velocity 

%% ED FIGURE 8D
plotModelingFigs(26);       %plot all fixed points in different regimes for velocities below threshold

%% ED FIGURE 8F
plotModelingFigs(27);       %plot all fixed points in different regimes for velocities above threshold


%---- INDEX KEY FOR PLOTTING MODELING RESULTS: PLOTMODELINGFIGS(INDEX) ---%
% FIGURE 2
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

%to generate all modeling figures, run: plotModelingFigs(1:27)
%-------------------------------------------------------------------------%
