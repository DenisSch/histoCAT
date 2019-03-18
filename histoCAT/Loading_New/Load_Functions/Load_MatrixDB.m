function [Sample_Set_arranged,Tiff_all,Tiff_name] = Load_MatrixDB(samplefolders,Sample_Set_arranged, ometif_name)
% LOAD_MATRIXDB: Main function for loading tiffs and mask data
%
% Input variable:
% samplefolders --> paths to the selected sample folders
% Sample_Set_arranged --> paths to all sample folders in session (historical)
% Mask_all --> segmentation masks of all samples (matrices)
%
% Output variables:
% Sample_Set_arranged --> paths to all sample folders in session (historical)
% Mask_all --> segmentation masks of all samples (matrices)
% Tiff_all --> tiff matrices of all samples (images / channels)
% Tiff_name --> tiff names of all samples (image / channel names)
%
% Histology Topography Cytometry Analysis Toolbox (histoCAT)
% Denis Schapiro - Bodenmiller Group - UZH

%If not stored before, add the sample folders to Sample_Set_arranged
[Sample_Set_arranged]= unique([Sample_Set_arranged samplefolders],'stable');

%Display update
disp('Samples arranged..Loading Masks,Tiffs of all Samples..');
hWaitbar = waitbar(0,'Loading all files from Database. This may take a while...');

%Loop through all the samples in DB
for i=1:size(Sample_Set_arranged,2)
    
    %Do not load mask and tiff, if tiff is multipage OME format (HARDCODED for
    %t-CycIF OME export for now)
    
    %Get all the files in the sample folder
    fileList = getAllFiles(char(Sample_Set_arranged(i)));
    
    %Extract tiffs (besides the one representing a mask)
    tiff_position = find(~cellfun('isempty',regexp(fileList,ometif_name)))';
    %tiff_position=1;
    %Extract image info from tiff
    Image_info = imfinfo(char(fileList(tiff_position)));
    
    %Load all mask from the TMA
    %[Mask_all] = Load_mask(Sample_Set_arranged,Mask_all,i);
    
    %Exclude from loading if multipage tiff
    if size(Image_info,1)>2        
       [Tiff_all,Tiff_name] = Load_multipage_tiff(Sample_Set_arranged,tiff_position,i);
    else
        %Load all tiff files and names for each image
        [Tiff_all,Tiff_name] = Load_tiff(Sample_Set_arranged,i);
    end
    
    %Update waitbar
    waitbar(i/size(Sample_Set_arranged,2), hWaitbar);
end

close(hWaitbar);

end

