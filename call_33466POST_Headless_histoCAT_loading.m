clear all
clc

%Specify actual parameters to pass to Headless_histoCAT_loading.m
%samplefolders_str = '/Users/ediknovik/Dropbox/Harvard-University/HMS-HiTS-LSP-Internship-2019/Triplet-Code/Melanoma_Triplet/auto-fluorescence-correction/';
samplefolders_str = '/home/en100/Images/';
tiff_name = '33466POST.ome.tif';
segmentationfolder_str = samplefolders_str;
%mask_name = '33466POSTcellMask.tif';
mask_name = 'cellMask.tif';
Marker_CSV = 'Triplet_40_markers.csv';
expansionpixels = '30';

tic
%Execute Headless_histoCAT_loading.m
Headless_histoCAT_loading...
    (samplefolders_str,...
    tiff_name,...
    segmentationfolder_str,...
    mask_name,...
    Marker_CSV,expansionpixels)
toc
disp('ran histoCAT!')

%Define path where marker mean .mat files are saved
mean_path = strcat(samplefolders_str,'mean_output/33466POST/');
CSV_33466POST= readtable(fullfile(mean_path, '33466POST.csv'));

%Extract marker names
Marker_list = table2array(readtable(Marker_CSV,'ReadVariableNames',false));
numMarkers = length(Marker_list);

%Save pixel arrays for each cell across markers
pixels_across_markers = {};

pixel_path = strcat(samplefolders_str,'pixel_output/33466POST/');

for i = 1:numMarkers
    marker_pixels = load(strcat(pixel_path,'Cell_33466POST',Marker_list{i,1},'.mat'),'get_pixels');
    
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

% % % save('pixels_across_markers.mat','33466_POST_pixels_across_markers','-v7.3');
% % % disp('saved patient pixels');

% % % numCells = size(pixels_across_markers{1,1},1);
% % % 
% % % %Perform pixel correlations for each cell across markers, pairwise
% % % corrs = zeros(numCells,3,3);
% % % pears_corrs = {}; %Save pixel correlations here
% % % 
% % % tic
% % % %Choose which marker channels to include for correlation matrices
% % % %AF channels:
% % % n_start = 2; 
% % % n_end = 4;
% % % 
% % % %Cycle 1 channels (excluding Hoechst):
% % % m_start = 6;
% % % m_end = 8;
% % % 
% % % for j = 1:numCells %for each cell
% % %     %two for loops for pairwise correlations across markers
% % %     %find correlations of AF markers and cycle 1
% % %     for n = n_start:n_end 
% % %         for m = m_start:m_end
% % %             
% % %             cell_j_marker_n = double(pixels_across_markers{1,n}{j,1});
% % %             cell_j_marker_m = double(pixels_across_markers{1,m}{j,1});
% % %             
% % %             pearson_corr = corrcoef(cell_j_marker_n, cell_j_marker_m);
% % %             %extract correlation coefficient
% % %             R = pearson_corr(1,2); %either off-diagonal (1,2) or (2,1) works
% % %             
% % %             %save correlation coefficient in matrix
% % %             corrs(j,n-1,m-5) = R;
% % %         end
% % %     end
% % %     %Reshape to 2D 3x3 matrix 
% % %     pears_corrs{j} = reshape(corrs(j,:,:),[3,3])';
% % % end
% % % toc
% % % 
% % % %Plot single cells with correlation overlay 
% % % corr_label = [];
% % % for k = 1:size(pears_corrs,2)
% % %     corr_label = vertcat(corr_label,diag(pears_corrs{k})');
% % % end
% % % 
% % % % figure
% % % % A488 = imread('Example.tif',2);
% % % % imagesc(A488)
% % % % figure
% % % % pERK = imread('Example.tif',6);
% % % % imagesc(pERK)
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter),exampleCSV.Y_position(~cells_to_filter),[],corr_label(:,1), 'filled')
% % % colorbar
% % % caxis([-1 1])
% % % title(strcat('Pixel Pearson Correlation Across Cells (A488-pERK)'))
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter),exampleCSV.Y_position(~cells_to_filter),[],corr_label(:,2), 'filled')
% % % colorbar
% % % caxis([-1 1])
% % % title(strcat('Pixel Pearson Correlation Across Cells (A555-AXL)'))
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter),exampleCSV.Y_position(~cells_to_filter),[],corr_label(:,3), 'filled')
% % % colorbar
% % % caxis([-1 1])
% % % title(strcat('Pixel Pearson Correlation Across Cells (A647-MITF)'))
% % % 
% % % %Plot pixel intensity of AF channels
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleA488(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'A488 intensity';
% % % caxis([7 11])
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleA555(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'A555 intensity';
% % % caxis([7 11])
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleA647(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'A647 intensity';
% % % caxis([7 11])
% % % 
% % % %Plot pixel intensity of 3 markers
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExamplepERK(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'pERK intensity';
% % % caxis([7 11])
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleAXL(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'AXL intensity';
% % % caxis([7 11])
% % % 
% % % figure
% % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleMITF(~cells_to_filter), 'filled')
% % % title('Example: Single cell position')
% % % xlabel('x position')
% % % ylabel('y position')
% % % c = colorbar;
% % % c.Label.String = 'MITF intensity';
% % % caxis([7 11])
% % % 
% % % % %Plot scatter plots with correlation overlay
% % % % figure
% % % % scatter(exampleCSV.X_position(~cells_to_filter), exampleCSV.Y_position(~cells_to_filter), [], exampleCSV.Cell_ExampleA488(~cells_to_filter), 'filled')
% % % % hold on;
% % % % scatter(exampleCSV.X_position(~cells_to_filter),exampleCSV.Y_position(~cells_to_filter),25,corr_label(:,1))
% % % % 
% % % 
% % % %Plot mean intensities for AF and Cycle 1 
% % % figure
% % % scatter(exampleCSV.Cell_ExamplepERK, exampleCSV.Cell_ExampleA488)
% % % xlabel('pERK')
% % % ylabel('A488')
% % % title('Mean Intensity')
% % % 
% % % figure
% % % scatter(exampleCSV.Cell_ExampleAXL, exampleCSV.Cell_ExampleA555)
% % % xlabel('AXL')
% % % ylabel('A555')
% % % title('Mean Intensity')
% % % 
% % % figure
% % % scatter(exampleCSV.Cell_ExampleMITF, exampleCSV.Cell_ExampleA647)
% % % xlabel('MITF')
% % % ylabel('A647')
% % % title('Mean Intensity')
% % % 
% % % %Plot single cell pixel scatter 
% % % cell = 10;
% % % marker1 = 8; %6-8
% % % marker2 = 4; %2-4
% % % 
% % % figure
% % % scatter(pixels_across_markers{1,marker1}{cell,1}, pixels_across_markers{1,marker2}{cell,1})
% % % xlabel(Marker_list(marker1))
% % % ylabel(Marker_list(marker2))
% % % title(strcat('Cell',num2str(cell),': Pixel Locations'))

%Plot scatter of all pixels for A488 and pERK

%Extract Hoechst1 and Hoechst2 marker pixels
Hoechst1_pixels = [];
Hoechst2_pixels = [];

for p = 1:numCells
    Hoechst1_pixels = vertcat(Hoechst1_pixels,pixels_across_markers{1,1}{p,1});
    Hoechst2_pixels = vertcat(Hoechst2_pixels,pixels_across_markers{1,5}{p,1});
end
Hoechst1_pixels = double(Hoechst1_pixels);
Hoechst2_pixels = double(Hoechst2_pixels);
%Pearson correlation
%corrcoef(Hoechst1_pixels,Hoechst2_pixels)

figure
scatter(log(Hoechst1_pixels),log(Hoechst2_pixels),5)
xlabel(strcat(Marker_list(1),{' '},'(log)'))
ylabel(strcat(Marker_list(5),{' '},'(log)'))
title(strcat('Pixel Scatter'))

%Extract A488 and pERK marker pixels
A488_pixels = [];
pERK_pixels = [];

for p = 1:numCells
    A488_pixels = vertcat(A488_pixels,pixels_across_markers{1,2}{p,1});
    pERK_pixels = vertcat(pERK_pixels,pixels_across_markers{1,6}{p,1});
end
A488_pixels = double(A488_pixels);
pERK_pixels = double(pERK_pixels);
%Pearson correlation
%corrcoef(A488_pixels,pERK_pixels)

figure
scatter(log(A488_pixels),log(pERK_pixels),5)
xlabel(strcat(Marker_list(2),{' '},'(log)'))
ylabel(strcat(Marker_list(6),{' '},'(log)'))
title(strcat('Pixel Scatter'))

%Extract A555 and AXL marker pixels
A555_pixels = [];
AXL_pixels = [];

for p = 1:numCells
    A555_pixels = vertcat(A555_pixels,pixels_across_markers{1,3}{p,1});
    AXL_pixels = vertcat(AXL_pixels,pixels_across_markers{1,7}{p,1});
end
A555_pixels = double(A555_pixels);
AXL_pixels = double(AXL_pixels);
%Pearson correlation
%corrcoef(A555_pixels,AXL_pixels)

figure
scatter(log(A555_pixels),log(AXL_pixels),5)
xlabel(strcat(Marker_list(3),{' '},'(log)'))
ylabel(strcat(Marker_list(7),{' '},'(log)'))
title(strcat('Pixel Scatter'))

%Extract A647 and MITF marker pixels
A647_pixels = [];
MITF_pixels = [];

for p = 1:numCells
    A647_pixels = vertcat(A647_pixels,pixels_across_markers{1,4}{p,1});
    MITF_pixels = vertcat(MITF_pixels,pixels_across_markers{1,8}{p,1});
end
A647_pixels = double(A647_pixels);
MITF_pixels = double(MITF_pixels);
%Pearson correlation
%corrcoef(A647_pixels,MITF_pixels)

figure
scatter(log(A647_pixels),log(MITF_pixels),5)
xlabel(strcat(Marker_list(4),{' '},'(log)'))
ylabel(strcat(Marker_list(8),{' '},'(log)'))
title(strcat('Pixel Scatter'))











