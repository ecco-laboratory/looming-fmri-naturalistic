%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS A SINGLE TASK
% aggregating over subject-level 1st-level SPM beta or contrast statmaps
% each of these should come in from targets calling script

% paths_nifti_betas: derived from target. 
% cell array of paths to the same beta/contrast nifti but different subjects
% paths_nifti_mean: optional, derived from target. 
% cell array of paths to the mean signal contrast nifti but different subjects
% if exists, will be divided out from the betas to yield percent signal change units
% region: char in canlabcore atlas syntax. I'm committing to the bit--this script again runs only over a single ROI.
% out_path: for the table of in-mask voxelwise betas by subject. one char array for one output file.
% recommended to put the ROI in the file names as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 

%% PREP ROI MASK
% use custom function that loads brainstemnavigator ROI for Bstem_SC
mask = select_this_atlas_subset(region);

%% GET SUBJECT NUMBERS FROM PATHS_NIFTI
% so it will be robust to potential differences in fold structure/order
subj_nums = zeros(length(paths_nifti_betas), 1);
for i=1:length(paths_nifti_betas)
    path_parts = strsplit(paths_nifti_betas{i}, "/");
    % THIS IS HARD CODED TO WORK ON SPM FIRST LEVEL BETAS SAVED IN OUR SPECIFIC FOLDER STRUCTURE
    subj_parts = strsplit(path_parts{10}, "-");
    subj_num = str2double(subj_parts{2});
    subj_nums(i) = subj_num;
end

%% LOAD STACK OF FIRST-LEVEL STATMAPS

contrasts = fmri_data(paths_nifti_betas);
contrasts_masked = apply_mask(contrasts, mask);
% transpose fmri_data.dat to get subject on the row and voxel on the column
contrasts_matrix = contrasts_masked.dat';

%% IF MEAN SIGNAL IS SPECIFIED, RESCALE STATMAPS TO PERCENT SIGNAL CHANGE
if exist('paths_nifti_mean', 'var') == 1
    meansignal = fmri_data(paths_nifti_mean);
    meansignal_masked = apply_mask(meansignal, mask);
    meansignal_matrix = meansignal_masked.dat';
    contrasts_matrix = contrasts_matrix ./ meansignal_matrix;
    contrasts_matrix = contrasts_matrix .* 100;
end

%% WRITE OUT MASKED BETAS TO TABLE
% write to tabular text data so you can read it into R for stats and graphs
writematrix([subj_nums, contrasts_matrix], out_path)
