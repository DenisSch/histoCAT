function [get_mean,get_mean_name] = Get_mean_batch(marker_position,...
    sessionData_name, sessionData_folder)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Load session
load(sessionData_name)

% Get global variables
global Sample_Set_arranged
global Mask_all

%Get current mask
Current_Mask = Mask_all.Image;

%get mean expression for multipage tiff
global ometif_name
large_tiff_location = fullfile(Sample_Set_arranged{1,1},ometif_name);
get_mean = struct2array(regionprops(Current_Mask, ...
    imread(large_tiff_location,marker_position), 'MeanIntensity'))';
%Check to make sure it is a string
get_mean_name_onlyname = (table2cell(Marker_list(marker_position,1)));
get_mean_name = strcat('Cell_',get_mean_name_onlyname(1,1));

save(fullfile(sessionData_folder,strcat(char(get_mean_name),'.mat')),...
    'get_mean','get_mean_name');

end

