function Workflow_DIFCO_cell_counting_cluster2_LaVision (nucleus_data_dir, Int_data_dir, threshold, area_dir, resolution, out_dir)

% Scripts for DIFCO project by Uhlen lab
% Ver 1.00
% These scripts are written by Shigeaki Kanatani 
%
% Contact: Per Uhlén, per.uhlen@ki.se
%          Shigeaki Kanatani, shigeaki.kanatani@ki.se

%% Summary
% nucleus_data is binary tif file directory which is going to be counted
% Int_data_dir is 16 bit tif file directory (usually in situ hybridization data)
% threshold is the value you will cut-off in Int_data_dir
% area_dir is the binary tif file directory which is segmented for whole sample shape.
% out_dir is the directory you will have data


%% Counting cells and basic statistic 
out_dir1 = [out_dir '\' 'Cell_Segment_Intensity'];
mkdir (out_dir1);
% Save data for counting data and summary

F_XYZID = Hmaxima2PointCloud_3D(nucleus_data_dir, out_dir1); 
% Count segmented cells (centroids) using bwconncomp and regionprops3

SE_shape= strel('disk', 8,0);
se=SE_shape.Neighborhood; 
%se is 2 dimension only. You can change if you want.
%This parameter is for the COLM system, 10X lens.

F_XYZIDint = Point3Dto2D_Dilate_MeanInt_Par (Int_data_dir, F_XYZID, se, out_dir1);
% Generate circle area defined with se (above) in each cell and measure mean intensity in Int_data_dir.

F_XYZIDint_Space = Resolution_adjustment_adjustable(F_XYZIDint, out_dir1, resolution);
% Adjustment of pixel to actual volume
% [0.585 0.585 5] is resolution specific to the COLM light sheet system 

[F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower] = Thresholding_PointCloud(F_XYZIDint_Space, threshold, out_dir1);
% This function outputs the point population which is more than threshold and less than threshold.
 
area_image = ImageSeries_openBinary(area_dir);
Volume3D = length(find(area_image == 1)) * resolution(1) * resolution(2) * resolution(3); 
%[0.585 0.585 5] is resolution specific to the COLM light sheet system 
%Volume of space used for analysis is calculated.

Basic_summary(F_XYZIDint_Space, F_XYZIDint_Space_ThUpper, F_XYZIDint_Space_ThLower, Volume3D, threshold, out_dir1);
% This save the summary in out_dir. 
% This script generate summary csv file of 
% ["Mean intensity"; "Median intensity"; "All cell number"; ...
% "Positive cell number"; "Negative cell number"; ...
% "positive cell ratio"; "total density of positive cell(mm3)"; ...
% "Volume(um3)"; "Threshold"];


%% Nearest nighbor distance analysis

% Nearest neighbor distance analysis
trial_number = 100; % The number of simulation you will do.
out_dir2 = [out_dir '\' 'Simulation']; % Save directory
mkdir (out_dir2);
mean_NN_MonteCarlo_3D_pdist_LaVision(area_dir, F_XYZIDint, threshold, trial_number, resolution, out_dir2)

 % This script perform MonteCarlo simulation for distribution of cells in the space defined by area_dir. 
 % The number of voxel of spaces are randomly selected and average nearest distance is obtained  
 % Threshold is used to extract points more than (>) thresold intensity (posivie cells).

 % Nearest neibor index(NNI) is to show how dense or sparse cells are distributed.
 % If NNI is less than 1, cells distribute densly more than simulation.
 % If NNI is more than 1, cells distribute sparsely more than simulation.
 % The real sample distribution value is compared with the value of the
 % simulation and check statistic difference using z-test.

%% Density analysis

% Thresholding
F_XYZIDint_Space_Thre = F_XYZIDint_Space(:,5) > threshold; 
positive_cell_list = F_XYZIDint_Space(F_XYZIDint_Space_Thre, 1:4);
  % Have a positive cell list which have average intensity value more than threshold

% Density analysis
radius = 100; % search radius to find cells
out_dir3 = [out_dir '\' 'Density_' num2str(radius) 'radius'];
mkdir (out_dir3); % save directory
[numberOfpoints] = pointcloud_density (positive_cell_list(:,1:3), radius, out_dir3); 
positive_cellDensity_list = cat(2, positive_cell_list, numberOfpoints);      % X, Y, Z, ID and density

% High density cell extraction
Density_th = 700; % The number of cells surrounding in radius for cut off
positive_highDensity_list = cat(2, positive_cell_list, numberOfpoints);     % X, Y, Z, ID and density, but before threshoulding
HighDensity_pointlist = numberOfpoints >= Density_th; % Extract cell which has neighbor more than Density_th  
positive_highDensity_list = positive_highDensity_list(HighDensity_pointlist, :);
save ([out_dir3 '\positive_highDensity_list-' num2str(Density_th) 'cell.mat'], 'positive_highDensity_list');       

 % Here we get positive cell list and its density.
 % Cells with high density (in publication 400 cells is used which is
 % average value of all samples) is extracted and saved.


%% Clustering analysis dist 20
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir4 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir4);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir4)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir10 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir10);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir10)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir11 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir11);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir11)

%% Clustering analysis 20(all cells)
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir5 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir5);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir5)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir12 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir12);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir12)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 20; % The distance (micro meter) to make cells clustered.
out_dir13 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir13);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir13)

%% Clustering analysis 10
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir6 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir6);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir6)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir14 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir14);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir14)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir15 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir15);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir15)

%% Clustering analysis 10(all cells)
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir7 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir7);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir7)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir16 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir16);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir16)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 10; % The distance (micro meter) to make cells clustered.
out_dir17 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir17);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir17)

%% Clustering analysis 5
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir8 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir8);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir8)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir18 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir18);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir18)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir19 = [out_dir '\' 'Point_cloud_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir19);
Hot_spot_analysis_pcsegdist_sample(positive_highDensity_list, distance_num, cell_num, Volume3D, out_dir19)

%% Clustering analysis 5(all cells)
% cell 5

cell_num = 5; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir9 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir9);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir9)

% cell 10

cell_num = 10; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir20 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir20);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir20)

% cell 20

cell_num = 20; % The number of cells in each cluster. The cluster which has less than cell_num won't be counted.
distance_num = 5; % The distance (micro meter) to make cells clustered.
out_dir20 = [out_dir '\' 'Point_cloud_all_' num2str(cell_num) 'cell_' num2str(distance_num) 'distance_pcsegdist'];
mkdir (out_dir20);
Hot_spot_analysis_pcsegdist_sample(positive_cellDensity_list, distance_num, cell_num, Volume3D, out_dir20)




end