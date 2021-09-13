function F_XYZIDint_Space = Resolution_adjustment_adjustable(F_XYZIDint, out_dir, resolution)

% Scripts for DIFCO project by Uhlen lab
% Ver 1.00
% These scripts are written by Shigeaki Kanatani 
%
% Contact: Per Uhlén, per.uhlen@ki.se
%          Shigeaki Kanatani, shigeaki.kanatani@ki.se

%% Summary 
 % Adjustment of pixel to actual volume
 % [0.585 0.585 5] is resolution of the COLM light sheet system
 % Resolution of LaVision system can be changed 
 % F_XYZIDint is arrary of X, Y, Z, ID, Intensity.
 
    F_XYZIDint_Space = F_XYZIDint;
    F_XYZIDint_Space(:,3) = F_XYZIDint_Space(:,3)* resolution(3);
    F_XYZIDint_Space(:,2) = F_XYZIDint_Space(:,2)* resolution(2);
    F_XYZIDint_Space(:,1) = F_XYZIDint_Space(:,1)* resolution(1);
    
    save([out_dir '\' 'F_XYZIDint_Space.mat'],'F_XYZIDint_Space');
    
    Header ={'X', 'Y', 'Z', 'ID', 'Intensity'};    
    T = array2table(F_XYZIDint_Space, 'VariableNames',Header);
    writetable(T, [out_dir '\' 'F_XYZIDint_Space.csv']);
    
    pointData = F_XYZIDint_Space(:,1:3);
    save([out_dir '\' 'pointData.mat'],'pointData');
    
end