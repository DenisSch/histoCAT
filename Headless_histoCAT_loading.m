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
sessionData_mean_folder = fullfile('mean_output',tiff_name_raw{1,1});
mkdir(sessionData_mean_folder);
sessionData_mean_name = fullfile('mean_output',tiff_name_raw{1,1},...
    strcat(tiff_name_raw{1,1},'_session.mat'));
save(sessionData_mean_name,'-v7.3');
sessionData_mean_name = fullfile('mean_output',tiff_name_raw{1,1},...
    strcat(tiff_name_raw{1,1},'_session.mat'));

sessionData_pixel_folder = fullfile('pixel_output',tiff_name_raw{1,1});
mkdir(sessionData_pixel_folder);
sessionData_pixel_name = fullfile('pixel_output',tiff_name_raw{1,1},...
    strcat(tiff_name_raw{1,1},'_session.mat'));
save(sessionData_pixel_name,'-v7.3');
sessionData_pixel_name = fullfile('pixel_output',tiff_name_raw{1,1},...
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
    sessionData_mean_folder = fullfile('mean_output',tiff_name_raw{1,1});
    mkdir(sessionData_mean_folder);
    sessionData_mean_name = fullfile('mean_output',tiff_name_raw{1,1},...
        strcat(tiff_name_raw{1,1},'_session.mat'));
    disp('saving session')
    save(sessionData_mean_name,'-v7.3');
    disp('session saved')
    
    sessionData_pixel_folder = fullfile('pixel_output',tiff_name_raw{1,1});
    mkdir(sessionData_pixel_folder);
    sessionData_pixel_name = fullfile('pixel_output',tiff_name_raw{1,1},...
        strcat(tiff_name_raw{1,1},'_session.mat'));
    disp('saving session')
    save(sessionData_pixel_name,'-v7.3');
    disp('session saved')
%end

%% Parfor loop or submit to cluster
% Check if all files already exist

% get mean expression for multipage tiff

parfor i=1:size(Marker_list,1)
%for i=1:size(Marker_list,1)
    % Run locally
    [get_mean,get_mean_name] = Get_mean_batch(i,sessionData_mean_name,tiff_name_raw);
    
    %     % Submit to system
    %     cluster_command = 'sbatch -p short -c 1 -t 1:00:00 --mem=8000 ';
    %     command_to_submit_change = strcat('--wrap="matlab -nodesktop -r \"/home/ds230/histoCAT/histoCAT/Loading_New/DataProcessing/Get_mean_batch.m(',num2str(i),')\""');
    %     systems_call{i} = strcat(cluster_command,command_to_submit_change);
    
end

parfor i=1:size(Marker_list,1)
%for i=1:size(Marker_list,1)
    % Run locally
    [get_pixels,get_pixels_name] = Get_pixels_batch(i,sessionData_pixel_name,tiff_name_raw);

end

delete(gcp);

% Combine get_mean's
get_mean_all = [];
get_mean_name_all = {};
disp('combine all means')

% Combine get_pixels
get_pixels_all = [];
get_pixels_name_all = {};
disp('combine all pixels')

for k=1:size(Marker_list,1)
    % load all Markers and create "get_pixel"
    Name_to_load = fullfile(sessionData_pixel_folder,...
        strcat('Cell_',tiff_name_raw{1,1},table2cell(Marker_list(k,1)),'.mat'));
    load(char(Name_to_load));
    % Create matrix with
    get_pixels_all = [get_pixels_all,get_pixels];
    get_pixels_name_all{1,k} = strcat('Cell_',tiff_name_raw{1,1},char(table2cell(Marker_list(k,1))));
    
    % load all Markers and create "get_mean"
    Name_to_load = fullfile(sessionData_mean_folder,...
        strcat('Cell_',tiff_name_raw{1,1},table2cell(Marker_list(k,1)),'.mat'));
    load(char(Name_to_load));
    % Create matrix with
    get_mean_all = [get_mean_all,get_mean];
    get_mean_name_all{1,k} = strcat('Cell_',tiff_name_raw{1,1},char(table2cell(Marker_list(k,1))));
end
disp('all means combined')
disp('all pixels combined')

%% Run spatial
%Run single cell processing
disp('run spatial')
[Fcs_Interest_all] = Process_SingleCell_Tiff_Mask_batch(Tiff_all,Tiff_name,...
    Mask_all,Fcs_Interest_all,HashID,get_mean_all,get_mean_name_all,...
    sessionData_mean_name,expansionpixels);
disp('save CSV')
writetable(Fcs_Interest_all{1,1},...
    fullfile(sessionData_mean_folder, strcat(tiff_name_raw{1,1},'.csv')));

%% Run AF correlation
disp('ran histoCAT!')

%Define path where marker mean .mat files are saved
% mean_path = strcat(samplefolders_str,'mean_output/33466POST/');
% CSV_33466POST= readtable(fullfile(mean_path, '33466POST.csv'));
%mean_path = '/home/en100/Pixel_Correlation_Results/mean_output/Example/';
%csv_file = readtable(fullfile(mean_path, 'Example.csv'));

csv_file = readtable(fullfile(sessionData_mean_folder, strcat(tiff_name_raw{1,1},'.csv')));

%Extract marker names
Marker_list = table2array(readtable(Marker_CSV,'ReadVariableNames',false));
numMarkers = length(Marker_list);

%Save pixel arrays for each cell across markers
pixels_across_markers = {};

% pixel_path = strcat(samplefolders_str,'pixel_output/33466POST/');
%pixel_path = '/home/en100/Pixel_Correlation_Results/pixel_output/Example/';

for i = 1:numMarkers
%     marker_pixels = load(strcat(pixel_path,'Cell_33466POST',Marker_list{i,1},'.mat'),'get_pixels');
    marker_pixels = load(fullfile(sessionData_pixel_folder,strcat('Cell_',tiff_name_raw{1,1},Marker_list{i,1},'.mat')),'get_pixels');
    %% Filter out cells that have 10 pixels or less
    marker_pixels_filtered = marker_pixels;
    
    pixel_lengths = [];
    for k = 1:length(marker_pixels_filtered.get_pixels)
        pixel_lengths(k) = length(marker_pixels_filtered.get_pixels{k,1});
    end
    pixel_lengths = pixel_lengths';
    cells_to_filter = pixel_lengths <= 50;
    
    marker_pixels_filtered.get_pixels(cells_to_filter) = [];
    %%
    pixels_across_markers{i} = marker_pixels_filtered.get_pixels;
end

save('pixels_across_markers.mat',strcat(tiff_name_raw{1,1},'_pixels_across_markers'),'-v7.3');
disp('saved patient pixels');

numCells = size(pixels_across_markers{1,1},1);

%Perform pixel correlations for each cell across markers, pairwise
corrs = zeros(numCells,4,4);
pears_corrs = {}; %Save pixel correlations here

tic
%Choose which marker channels to include for correlation matrices
%AF channels:
n_start = 1; 
n_end = 4;

%Cycle 1 channels (including Hoechst):
m_start = 5;
m_end = 8;

for j = 1:numCells %for each cell
    %two for loops for pairwise correlations across markers
    %find correlations of AF markers and cycle 1
    for n = n_start:n_end 
        for m = m_start:m_end
            
            cell_j_marker_n = double(pixels_across_markers{1,n}{j,1});
            cell_j_marker_m = double(pixels_across_markers{1,m}{j,1});
            
            pearson_corr = corrcoef(cell_j_marker_n, cell_j_marker_m);
            %extract correlation coefficient
            R = pearson_corr(1,2); %either off-diagonal (1,2) or (2,1) works
            
            %save correlation coefficient in matrix
            corrs(j,n-n_start+1,m-m_start+1) = R;
        end
    end
    %Reshape to 2D 3x3 matrix
    pears_corrs{j} = reshape(corrs(j,:,:),[4,4]);
end
toc

%Plot single cells with correlation overlay 
corr_label = [];
for k = 1:size(pears_corrs,2)
    corr_label = vertcat(corr_label,diag(pears_corrs{k})');
end

tiff_name_raw = strsplit(tiff_name,'.');

f1 = figure('visible', 'off');
scatter(csv_file.X_position(~cells_to_filter),csv_file.Y_position(~cells_to_filter),[],corr_label(:,1), 'filled')
colorbar
caxis([-1 1])
title(strcat('Pixel Pearson Correlation Across Cells (Hoechst1-Hoechst2)'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-Hoechst1-Hoechst2-correlation','.fig'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-Hoechst1-Hoechst2-correlation','.tif')) 
close(f1)

f2 = figure('visible', 'off');
scatter(csv_file.X_position(~cells_to_filter),csv_file.Y_position(~cells_to_filter),[],corr_label(:,2), 'filled')
colorbar
caxis([-1 1])
title(strcat('Pixel Pearson Correlation Across Cells (A488-pERK)'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A488-pERK-correlation','.fig'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A488-pERK-correlation','.tif')) 
close(f2)

f3 = figure('visible', 'off');
scatter(csv_file.X_position(~cells_to_filter),csv_file.Y_position(~cells_to_filter),[],corr_label(:,3), 'filled')
colorbar
caxis([-1 1])
title(strcat('Pixel Pearson Correlation Across Cells (A555-AXL)'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A555-AXL-correlation','.fig'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A555-AXL-correlation','.tif')) 
close(f3)

f4 = figure('visible', 'off');
scatter(csv_file.X_position(~cells_to_filter),csv_file.Y_position(~cells_to_filter),[],corr_label(:,4), 'filled')
colorbar
caxis([-1 1])
title(strcat('Pixel Pearson Correlation Across Cells (A647-MITF)'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A647-MITF-correlation','.fig'))
saveas(gcf,strcat(tiff_name_raw{1,1},'-A647-MITF-correlation','.tif')) 
close(f4)


end

