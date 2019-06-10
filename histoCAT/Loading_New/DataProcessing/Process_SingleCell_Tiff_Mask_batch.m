function [Fcs_Interest_all] = Process_SingleCell_Tiff_Mask_batch(...
    Tiff_all,Tiff_name,Mask_all,Fcs_Interest_all,HashID,...
    get_mean_all,get_mean_name_all,sessionData_name,expansionpixels,varargin)
% PROCESS_SINGLECELL_TIFF_MASK_BATCH:
% This function stores the single cell information along with the expansion
% of pixels in mask mentioned by user and gets their neighbours. All
% masks/tiffs will be used to store these single cell information
% into the matrix of tables 'Fcs_Interest_all' of each ImageId/HashID processed.
%
% This is optimized for batch processing
%
% Input variables:
% Mask_all --> segmentation masks of all samples
% Tiff_all --> tiff matrices of all samples (images / channels)
% Tiff_name --> tiff names of all samples (image / channel names)
% HashID --> Unique folder IDs (!!!GLOBAL!!!)
% expansionfeature --> amount of pixels to expand from each cell in order
% to look for neighboring cells (set by user)
%
% Output variables;
% Fcs_Interest_all --> data in fcs format (first column: ImageID, second
% column: CellID, third column: marker1, etc.)
%
% Histology Topography Cytometry Analysis Toolbox (histoCAT)
% Denis Schapiro - Bodenmiller Group - UZH

% Load session
load(sessionData_name)

% Get global variables
global samplefolders

%If session was loaded for the first time, ask for the pixelexpansion to
%use in order to search for neighboring cells

% expansion = {'4','4'};
% put('expansionfeature',expansion(1));
% put('expansion_range',str2num(cell2mat(expansion(2))));

%Initialize variables
length_neighbr = [];
sizes_neighbrs = [];
idx_cel = find(~cellfun('isempty',struct2cell(Mask_all)));
allvarnames_nospatial = get_mean_name_all;

%Add spatial features to variable names
BasicFeatures = {'Area', 'Eccentricity', 'Solidity', 'Extent', ...
    'EulerNumber', 'Perimeter',...
    'MajorAxisLength', 'MinorAxisLength', 'Orientation'};
%Add X and Y
XY = {'X_position','Y_position'};
allvarnames = [allvarnames_nospatial, BasicFeatures,XY];

%Get current mask
Current_Mask = Mask_all.Image;

%get mean expression for multipage tiff
global tiff_name

% Load all markers

large_tiff_location = fullfile(samplefolders{1,1},tiff_name);
Current_singlecellinfo_nospatial = log(get_mean_all);

tic
%Get spatial features similar to CellProfiler
props_spatial = regionprops(Current_Mask, BasicFeatures(~strcmp(BasicFeatures,'FormFactor')));
%Add X and Y coordinates to output
props_spatial_XY = regionprops(Current_Mask, 'Centroid');
toc

%Add spatial information to data matrix: variable names and
%data
global Marker_list
Current_channels_nospatial = Marker_list{:,1}';
Current_channels = [Current_channels_nospatial, BasicFeatures,XY];

BasicFeatures_Matrix = [cat(1,props_spatial.Area),...
    cat(1,props_spatial.Eccentricity),...
    cat(1,props_spatial.Solidity),...
    cat(1,props_spatial.Extent),...
    cat(1,props_spatial.EulerNumber),...
    cat(1,props_spatial.Perimeter),...
    cat(1,props_spatial.MajorAxisLength),...
    cat(1,props_spatial.MinorAxisLength),...
    cat(1,props_spatial.Orientation),...
    cat(1,props_spatial_XY.Centroid)];

save('Save_before_bug.mat','Current_singlecellinfo_nospatial','BasicFeatures_Matrix','-v7.3');
disp('saved');

Current_singlecellinfo= [Current_singlecellinfo_nospatial, BasicFeatures_Matrix];

%Function call to expand cells and get the neighbrcellIds
[ Fcs_Interest_all,length_neighbr,sizes_neighbrs ] = NeighbrCells_histoCATsinglecells(1,allvarnames,Current_channels,Current_Mask,Current_singlecellinfo,...
    Fcs_Interest_all,length_neighbr,sizes_neighbrs,HashID,expansionpixels);

end
