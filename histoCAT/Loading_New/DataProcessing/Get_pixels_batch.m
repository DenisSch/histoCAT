function [get_pixels,get_pixels_name] = Get_pixels_batch(marker_position,...
    sessionData_name,tiff_name_raw)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Load session
load(sessionData_name)

clear get_mean

% Get global variables
global samplefolders
global Mask_all

%Get current mask
Current_Mask = Mask_all.Image;

% get mean expression for multipage tiff
global tiff_name
large_tiff_location = fullfile(samplefolders{1,1},tiff_name);
    %get_pixels = struct2array(regionprops(Current_Mask, ...
    get_pixels = struct2cell(regionprops(Current_Mask, ...    
        imread(large_tiff_location,marker_position), 'PixelValues'))';
%Check to make sure it is a string
get_pixels_name_onlyname = (table2cell(Marker_list(marker_position,1)));
get_pixels_name = strcat('Cell_',tiff_name_raw{1,1},get_pixels_name_onlyname(1,1));

save(fullfile(sessionData_folder,strcat(char(get_pixels_name),'.mat')),...
    'get_pixels','get_pixels_name');

end

