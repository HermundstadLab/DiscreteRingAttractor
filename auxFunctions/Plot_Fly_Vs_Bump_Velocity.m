%% A function for plotting the linear fit between the fly and bump velocities

function Plot_Fly_Vs_Bump_Velocity(FlyAV_Good, Bump_AV_Good, LinearFit_Right, LinearFit_Left, TrialInd, PlotDir)
    
figure; 
hold on; 
scatter(FlyAV_Good, Bump_AV_Good, 10, 'MarkerFaceColor', 'k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha', .02, 'MarkerEdgeAlpha', .02)
plot([0:1:350],  [0:1:350]*LinearFit_Right(1) + LinearFit_Right(2),'r')
plot([-350:1:0], [-350:1:0]*LinearFit_Left(1) + LinearFit_Left(2),'b')
xlim([-300 300])
ylim([-300 300])
xlabel('Fly AV (deg/s)')
ylabel('Bump AV (deg/s)')
title(['Fly ' num2str(TrialInd)])

set(gcf,'PaperOrientation','portrait','PaperPosition',[0 0 6 5.5],'PaperSize',[6 5.5]) % Vectors are width then height
print(gcf, '-dpng', '-r350', [PlotDir '\Fly_' num2str(TrialInd) '_VelocityGain'])
close gcf




