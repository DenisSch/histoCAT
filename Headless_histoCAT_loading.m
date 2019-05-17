function [] = Headless_histoCAT_loading...
    (samplefolders_str,...
    tiff_name,...
    segmentationfolder_str,...
    mask_name,...
    Marker_CSV,expansionpixels)
%HEADLESS_HISTOCAT_LOADING Headless loading for histoCAT
%   This function enables headless loading in histoCAT. This also enables
%   O2 cluster processing. It is optimized for large multipage tiffs.

% Histology Topography Cytometry Analysis Toolbox (histoCAT)
% Denis Schapiro - Independent Fellow -  Harvard and Broad Institute - 2019
addpath(genpath(pwd))


%% Please adapt this part to your data
%Load multipage tiff file(s)
%samplefolders_str = '/Users/denis/Desktop/Test_folder';
samplefolders = {samplefolders_str};
%tiff_name = '33466POST.ome.tif';
tiff_name_raw = strsplit(tiff_name,'.');

% Load mask
%segmentationfolder_str = '/Users/denis/Desktop/Test_folder'
%mask_name = '33466POST_cellMask.tif';
mask_location = fullfile(segmentationfolder_str,mask_name);

% Where is the marker list
%Marker_CSV = '/Users/denis/Desktop/Test_folder/Triplet_40_markers.csv';
Marker_list = readtable(Marker_CSV,'ReadVariableNames',false);

% global just for batch mode
global Marker_list
global transform_option_batch
global Tiff_name

% Define pixel expansion
% expansionpixels = 30;

% Transformation: option_list = {'Do not transform data','arcsinh','log'};
transform_option_batch = 'log';

%% Check if data already exist otherwise extract extract code from "Master_LoadSamples"
% Call global variables
global Mask_all
global Fcs_Interest_all
global HashID

% Function call to store the sample folder
[samplefolders,fcsfiles_path,HashID] = Load_SampleFolders(HashID,samplefolders);

% Load all the db files
[Mask_all,Tiff_all,...
    Tiff_name]= Load_MatrixDB_batch(mask_location);

% Save session to folder
sessionData_folder = fullfile('output',tiff_name_raw{1,1});
mkdir(sessionData_folder);
sessionData_name = fullfile('output',tiff_name_raw{1,1},...
    strcat(tiff_name_raw{1,1},'_session.mat'));
save(sessionData_name,'-v7.3');
sessionData_name = fullfile('output',tiff_name_raw{1,1},...
    strcat(tiff_name_raw{1,1},'_session.mat'));

% if isfile(sessionData_name)
%     disp('sessionData found');
%     load(sessionData_name);
%     disp('sessionData loaded');
% else
    %Call global variables
    global Mask_all
    global Fcs_Interest_all
    global HashID
    
    % Function call to store the sample folder
    [samplefolders,fcsfiles_path,HashID] = Load_SampleFolders(HashID,samplefolders);
    
    % Load all the db files
    [Mask_all,Tiff_all,...
        Tiff_name]= Load_MatrixDB_batch(mask_location);
    
    % Save session to folder
    sessionData_folder = fullfile('output',tiff_name_raw{1,1});
    mkdir(sessionData_folder);
    sessionData_name = fullfile('output',tiff_name_raw{1,1},...
        strcat(tiff_name_raw{1,1},'_session.mat'));
    disp('saving session')
    save(sessionData_name,'-v7.3');
    disp('session saved')
%end

%% Parfor loop or submit to cluster
% Check if all files already exist

% get mean expression for multipage tiff

parfor i=1:size(Marker_list,1)
    % Run locally
    [get_mean,get_mean_name] = Get_mean_batch(i,sessionData_name,tiff_name_raw);
    
    %     % Submit to system
    %     cluster_command = 'sbatch -p short -c 1 -t 1:00:00 --mem=8000 ';
    %     command_to_submit_change = strcat('--wrap="matlab -nodesktop -r \"/home/ds230/histoCAT/histoCAT/Loading_New/DataProcessing/Get_mean_batch.m(',num2str(i),')\""');
    %     systems_call{i} = strcat(cluster_command,command_to_submit_change);
    
end

delete(gcp);

% Combine get_mean's
get_mean_all = [];
get_mean_name_all = {};
disp('combine all means')
for k=1:size(Marker_list,1)
    % load all Markers and create "get_mean"
    Name_to_load = fullfile(sessionData_folder,...
        strcat('Cell_',tiff_name_raw{1,1},table2cell(Marker_list(k,1)),'.mat'));
    load(char(Name_to_load));
    % Create matrix with
    get_mean_all = [get_mean_all,get_mean];
    get_mean_name_all{1,k} = strcat('Cell_',tiff_name_raw{1,1},char(table2cell(Marker_list(k,1))));
end
disp('all means combined')

%% Run spatial
%Run single cell processing
disp('run spatial')
[Fcs_Interest_all] = Process_SingleCell_Tiff_Mask_batch(Tiff_all,Tiff_name,...
    Mask_all,Fcs_Interest_all,HashID,get_mean_all,get_mean_name_all,sessionData_name);
disp('save CSV')
writetable(Fcs_Interest_all{1,1},...
    fullfile(sessionData_folder, strcat(tiff_name_raw{1,1},'.csv')));

end

