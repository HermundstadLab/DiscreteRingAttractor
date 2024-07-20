% This function bins the bump velocity residuals by the bump's position

function [Binned] = Bin_Bump_Velocity_Residuals(AllFly, HD_Bin_Edges)


% Loop over flies and do the binning
Binned.Left_N=[];   Binned.Bump_HD_Good_Left=[];   Binned.Bump_AV_Good_Residual_Left=[];   Binned.Bump_AV_Good_Residual_Left_STD=[];
Binned.Right_N=[];  Binned.Bump_HD_Good_Right=[];  Binned.Bump_AV_Good_Residual_Right=[];  Binned.Bump_AV_Good_Residual_Right_STD=[];
for fff = 1 : length(AllFly.Bump_HD_Good_Left)
    for bbb = 1 : length(HD_Bin_Edges)-1
        
        % Bin left turn data 
        left_logica =  AllFly.Bump_HD_Good_Left{fff} >=  HD_Bin_Edges(bbb) & AllFly.Bump_HD_Good_Left{fff} < HD_Bin_Edges(bbb+1); 
        Binned.Left_N(fff, bbb)                         =   sum(left_logica); 
        Binned.Bump_HD_Good_Left(fff,bbb)               =   mean(AllFly.Bump_HD_Good_Left{fff}(left_logica));
        Binned.Bump_AV_Good_Residual_Left(fff,bbb)      =   mean(AllFly.Bump_AV_Good_Residual_Left{fff}(left_logica));
        Binned.Bump_AV_Good_Residual_Left_STD(fff,bbb)  =   std(AllFly.Bump_AV_Good_Residual_Left{fff}(left_logica));

        % Bin right turn data 
        right_logica =  AllFly.Bump_HD_Good_Right{fff} >=  HD_Bin_Edges(bbb) & AllFly.Bump_HD_Good_Right{fff} < HD_Bin_Edges(bbb+1); 
        Binned.Right_N(fff, bbb)                        =   sum(right_logica); 
        Binned.Bump_HD_Good_Right(fff,bbb)              =   mean(AllFly.Bump_HD_Good_Right{fff}(right_logica));
        Binned.Bump_AV_Good_Residual_Right(fff,bbb)     =   mean(AllFly.Bump_AV_Good_Residual_Right{fff}(right_logica));
        Binned.Bump_AV_Good_Residual_Right_STD(fff,bbb) =   std(AllFly.Bump_AV_Good_Residual_Right{fff}(right_logica));

    end
end