%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% each of these should come in from targets calling script

% path_fmri_data: path to CSV input file with subjects (or whatever) on the rows and fmri_data whole brain voxel on the columns
% must have come from canlabtools cause it's going back into canlabtools
% path_fmri_data_2: optional, like path_fmri_data. if specified, will be loaded in and combined with the fmri_data in path_fmri_data 
% using the specified arithmetic operator (to calculate a voxelwise comparison before averaging voxels within parcel)
% fun_compare: if path_fmri_data_2 is specified, this should be the name of the comparator function as a string. goes into feval()
% out_path: full path to write CSV output file of parcel avgs.
% nifti name (e.g., subject) on the rows and parcel on the columns

%% load in fmri_data-compatible matrix and slap it onto a sample object
% set constant sample nifti using one subject's first-level beta image
% should be fine because the non-dynamic scan parameters of controlled and naturalistic are the same
% the logic below only works (?) when we patch the data into an image that was originally 3D, 
% so we can't do this with the first vol of a timeseries nifti
data = fmri_data('/home/data/eccolab/SPLaT_fMRI/ignore/models/task-naturalistic/acq-mb8/sub-0001/model-boxcar/smoothed-4mm/beta_0001.nii');

wb_data = readmatrix(path_fmri_data);

if exist('path_fmri_data_2', 'var') == 1
    wb_data_2 = readmatrix(path_fmri_data_2);
    wb_data = feval(fun_compare, wb_data, wb_data_2);
end

% transpose again, remember, fmri_data.dat needs voxels on the row
data.dat = wb_data';

%% loop over ROIs and average voxels within each
% select_atlas_subset expects a cell array
% can have more than one ROI and it will be the union of all the ROIs

% for a whole-brain-esque analysis, ROIs are now hard-set in this script!
% atm, go through all glasser parcels plus amygdala, superior colliculus, and thalamus pulvinar & LGN
% TODO: annoyingly, the canlab2018 atlas doesn't breakdown amygdala by hemisphere. 
% only by nuclei but each nucleus ROI has both hemispheres?
% newer canlab atlases do break it up by hemisphere but no longer by nucleus (prob fine)
% however, I can't figure out how to get canlab2023 to work with our existing download of Bianciardi Brainstem Navigator
atlas = load_atlas('canlab2018');
% bear in mind that this separates each of the Ctx parcels by hemisphere
rois = atlas.labels(startsWith(atlas.labels, 'Ctx'));
% but these are bilateral (SC has the option to split by hemi but by putting only the prefix canlabtools will group)
rois = [rois {'Amygdala'} {'Bstem_SC'} {'Thal_Pulv'} {'Thal_LGN'}];

% preallocate output array to have subject on the rows and parcel on the columns
parcel_means = zeros(height(wb_data), length(rois));

for i=1:length(rois)
    disp(['Current ROI: ' rois{i}])
    % critically, when it applies a mask to a stack of niftis
    % it will keep any voxel that is within the mask for ANY image
    % so we don't have to handle when some subjects have zero voxels in the mask
    masked_dat = apply_mask(data, select_this_atlas_subset(rois(i), atlas));
    % AVERAGE OVER VOXELS (ROW) TO GET A SINGLE VALUE FOR EACH PARCEL
    % then transpose to get file (subject, implicit) on the row and parcel on the column
    parcel_means(:, i) = mean(masked_dat.dat, 1)';
    
end

% apply the nifti paths as row names
parcel_means_table = array2table(parcel_means, VariableNames=rois);

%% SAVE OUT
% the parcel name will be stored in the file name
writetable(parcel_means_table, out_path);
