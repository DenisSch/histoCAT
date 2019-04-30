function [Mask_all,Tiff_all,Tiff_name] = Load_MatrixDB_batch(mask_location)
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

% Load mask with single precision 32 bit float
Mask_single_precision = imread(mask_location);
% Extract unique values corresponding to CellIDs
[a b c] = unique(Mask_single_precision);
% HARDCODED!!! First "CellID" is background zeros
a = [0:size(a,1)];
% Create a 32ubit mask without float
a_uint32 = uint32(a*(2^32));
% Reshape vector to matrix with unique value location
Mask_all(1).Image = reshape(a_uint32(c),size(Mask_single_precision));

% Process
[Tiff_all,Tiff_name] = Load_multipage_tiff(1);

end

