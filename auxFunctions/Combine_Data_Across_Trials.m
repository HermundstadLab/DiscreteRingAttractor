% This function combines data across trials for each fly

function [FlyData]=Combine_Data_Across_Trials(DarkData)

Flies=unique(DarkData.Fly);
FlyData.FlyForward_w2p=[];
FlyData.FlyRotate_w2p=[];
FlyData.BumpPosition=[];
FlyData.BumpR=[];
FlyData.FlyPosition_w2p_aligned=[];
FlyData.TwoP_Times=[];
FlyData.TwoP_Fs=[];
for fff=1:length(Flies)
    FlyLogica = find(DarkData.Fly == Flies(fff));
    FlyData.FlyForward_w2p(fff,:)           =  reshape( DarkData.FlyForward_w2p(FlyLogica,:)', 1, []);
    FlyData.FlyRotate_w2p(fff,:)            =  reshape( DarkData.FlyRotate_w2p(FlyLogica,:)', 1, []);
    FlyData.BumpPosition(fff,:)             =  reshape( DarkData.BumpPosition(FlyLogica,:)', 1, []);
    FlyData.BumpR(fff,:)                    =  reshape( DarkData.BumpR(FlyLogica,:)', 1, []);
    FlyData.FlyPosition_w2p_aligned(fff,:)  =  reshape( DarkData.FlyPosition_w2p_aligned(FlyLogica,:)', 1, []);
    FlyData.TwoP_Times(fff,:)               =  DarkData.TwoP_Times(FlyLogica(1),:);
    FlyData.TwoP_Fs(fff,:)                  =  DarkData.TwoP_Fs(FlyLogica(1),:); 
end


