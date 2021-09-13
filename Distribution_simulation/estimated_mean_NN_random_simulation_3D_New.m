function estimated_mean_NN_random_simulation_3D_New(area_dir, points, resolution, cell_size_diameter, trial_number, out_dir)

% script is improved version of nearest neighborhood analysis from DIIFCO
% paper in Nature Biomedical Engineering (Tanaka et al., 2020)

% This script considers the size of the simulated cell and all random cells
% has distance each other to reflect real situation. All cell volume is
% designed to be inside of the sample shape.

%% Save data
% data_ind_dists: Nearest neighborhood analysis of real segmented point analysis (points)
% ind_dists_All: Nearest neighborhood analysis of simulated data from area_dir
% emNN_mean: Mean values of nearest neighborhood distance in each simulation
% emNN_random_std: Standard deviation of nearest neighborhood distance from each simulation
% emNN_random_mean: Mean value of nearest neighborhood distance of all simulation
% points_space_record: All xyz coordinate of simulated points.
% record_dists_all: Selection history of simulated points by criteria of overlap and outside of area.

%% Contact
% Shigeaki Kanatani, Staff Scientist, Karolinska Institutet
% shigeaki.kanatani@ki.se

%% START
%--------------------------------------------------------------------------------
% area_dir is the folder path for binary area data, which has isotropic resolution
% points is xyz coordinate (3 columns) in real space (um)
% resolution should be isotropic and one value; In Csaba's case 5um isotropic / voxel; example 5
% cell_size_diameter is the sphere diameter and the value is a diameter of shpere(maicro-meter); example, 50
% trial_number is how many time you want to simulate
% out_dir is the folder path to save output data


% load area_image
area_image = ImageSeries_openBinary(area_dir);

% get area_rem point number
indM = find(area_image == 1);
area_point_number = length(indM);

% get the number of segmentated number
points = points(:,1:3);   % points should be the csv file of 3 columns, but just in case
point_number = length(points);

sz=size(area_image); %Image_size=[y,x,z]; %pixel value

% preallocation for recording (everything is zero or empty)
emNN_mean = zeros(trial_number,1);
points_space_record = cell(trial_number,1);
record_dists_all = cell(trial_number,1);
ind_dists_All = zeros(point_number,2, trial_number);

% generate perimeter (Area - 1 voxel erosion)
se = strel('sphere', 1);   % this make a structure element of sphere of diameter 3.
area_image_rem = imerode(area_image, se);   % this remove one voxel from the original sample shape
area_image_peri = area_image - area_image_rem;  % this process leave voxel of perimeter and remove inside.
indM_peri = find(area_image_peri == 1);         % this find the position of shape in the single indexing
area_peri_point_number = length(indM_peri);
[y_peri, x_peri, z_peri] = ind2sub(sz, indM_peri);             % convert peri_points into x, y, z coordinate
points_space_peri = [resolution * x_peri resolution * y_peri resolution * z_peri]; % Get points in real space
ptCloud_peri = pointCloud(points_space_peri);       % convert perimeter points into pointCloud

tic

poolobj = gcp('nocreate');
delete(poolobj);
no_of_workers = 12;
parpool ('local',no_of_workers);  

parfor t=1:trial_number
        
        % Preallocation
        ind_dists = zeros(point_number, 2);
        record_dists = [];
        points_space_re = [];
        points_space = [];
        point_logical = false(point_number,1); 
 
        % Precondition
        points_needs = point_number - sum(point_logical(:)) ;  % this should be the same value with point_number
        
        while points_needs > 0    % until point all suffice the condition
        
        % Getting random points from space
        random_points = datasample(indM, points_needs);     % Get random points from segmented area
        [y, x, z] = ind2sub(sz, random_points);             % convert random_points into x, y, z
        points_space_re = [resolution * x resolution * y resolution * z]; % Get points in real space
        
        points_space = cat(1, points_space(point_logical,:), points_space_re);  % replace bad points to new random points
        ptCloud = pointCloud(points_space);     % conversion all points to pointcloud              
        point_logical = false(point_number,1);                  % preallocation

        % remove points which area is overlapping or outside of the space
        for o=1:point_number
            % Calculate the distance to second nearest point in the random points 
            [indices_seg,dists_seg] = findNearestNeighbors(ptCloud, ptCloud.Location(o,:), 2);
            % Calculate the distance to nearest point of peri-meter point of surface
            [indices_peri,dists_peri] = findNearestNeighbors(ptCloud_peri, ptCloud.Location(o,:), 1);
            % if two point are closer than 50um, or a point is close more than 25um(half of cell size), remove it
            if dists_seg(2) < cell_size_diameter || dists_peri < cell_size_diameter/2  
    
                point_logical(o) = 0;       %if it is overlapped or very close to surface, value will be 0
            else
                point_logical(o) = 1;       %Otherwiese, it is 1(True)
            end 
                
            ind_dists(o,1) = double(dists_seg(2));         % For recording
            ind_dists(o,2) = double(dists_peri);
            
        end    
              
        points_needs = point_number - sum(point_logical(:));    % Calculate how many points need to be reobtained.
        
        record_dists = cat(2, record_dists, ind_dists, point_logical);
        % For recording. Order is cell-cell distance, cell-surface distance and flag(0 is not ok) 
        
        end
        
%% Calcurate nearest neighbor distance        
        for p=1:point_number
            ptCloud = pointCloud(points_space);
             % Calculate distance to second nearest point because first point is itself. 
            [indices, dists] = findNearestNeighbors(ptCloud, ptCloud.Location(p,:),2);  
            ind_dists(p,:) = [double(dists(2)), double(indices(2))]; % Take info about second points
            %disp(['Processing point=' num2str(p)]);
        
        end
            
        record_dists_all{t} = record_dists;             %for recording
        points_space_record{t} = points_space;          %for recodring
        emNN_mean(t) = mean(ind_dists(:,1));            %for recording
        ind_dists_All(:,:,t) = ind_dists;               %for recording
        disp(['Simulation trial=' num2str(t) ' done.']);
end

%% Calculation of the mean NN value of actual data
        data_ind_dists = zeros(point_number,2);
        y = points(:,1);
        x = points(:,2);
        z = points(:,3);
        data_points_space = [x y z];   % the points already um scale in real space
     
        parfor p=1:point_number
            
            ptCloud_seg = pointCloud(data_points_space);
            [indices, dists] = findNearestNeighbors(ptCloud_seg, ptCloud_seg.Location(p,:),2);
            data_ind_dists(p,:) = [double(indices(2)), double(dists(2))];
            disp(['Processing point=' num2str(p)]);
            
        end
        
        data_mean = mean(data_ind_dists(:,2));
        
%% Calculation of NN index
    
    emNN_random_mean = mean(emNN_mean);
    emNN_random_std = std(emNN_mean, 1);
    %enNN_random_stE = emNN_random_std / sqrt(point_number);

    [h,p,ci,zval] = ztest(data_ind_dists(:,1), emNN_random_mean, emNN_random_std);
    
    NN_index = data_mean/ emNN_random_mean;

        save ([out_dir '\' 'data_ind_dists.mat'],'data_ind_dists');       
        NNA_summary = cat(2, data_mean, emNN_random_mean, emNN_random_std, ... 
                                  point_number, trial_number, NN_index, ...
                                  h, p, ci(1), ci(2), zval);
        Header ={'Data_Mean_distance', 'Simulation_Mean_distance', 'Simulation_Standard_deviation', ...
                 'Point_number', 'Trial_number', 'NN_index', ...
                 'Hypothesis', 'p-value', 'Confidence interval1', 'Confidence interval2','z-value'};
        csvwrite_with_headers([out_dir '\' 'NNA_summary.csv'], NNA_summary,Header);
        save([out_dir '\' 'NNA_summary.mat'],'NNA_summary');

        save ([out_dir '\' 'ind_dists_All.mat'],'ind_dists_All');
        save ([out_dir '\' 'emNN_mean.mat'],'emNN_mean');
        save ([out_dir '\' 'emNN_random_mean.mat'],'emNN_random_mean');
        save ([out_dir '\' 'emNN_random_std.mat'],'emNN_random_std');

        save ([out_dir '\' 'points_space_record.mat'],'points_space_record');    
        save ([out_dir '\' 'record_dists_all.mat'],'record_dists_all');    

        
    mkdir ([out_dir '\' 'simulation_csv']);    
    for s=1:trial_number
        
        simulation_save = points_space_record{s};
        writematrix(simulation_save, [out_dir '\simulation_csv\simulation_csv_' num2str(s,'%.4d') '.csv']);
    end
        
poolobj = gcp('nocreate');
delete(poolobj);

toc

end
