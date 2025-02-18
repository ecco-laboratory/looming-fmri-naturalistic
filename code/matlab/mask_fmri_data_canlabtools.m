%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS A SINGLE SUBJECT/TASK
% with multiple runs per subject
% the nested cell arrays are one level deep only, by run
% each of these should come in from targets calling script

% paths_nifti: derived from target. 
% nested cell array of paths to smoothed timeseries niftis
% trs_to_use: vector of valid TR indices for SPM indexing into the 4D volume
% tr_duration: in seconds, for the band-pass filter preproc
% paths_confounds: derived from target. 
% nested cell array to txts of SPM-compatible 2D matrices of confound regressors of interest
% already had the initial volumes of each removed
% region: char in canlabcore atlas syntax. I'm committing to the bit--this script again runs only over a single ROI.
% out_path: for masked, preprocessed canlabtools-able multivoxel timecourse. one char array for one output file.
% recommended to put the ROI in the file names as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 

n_runs = length(paths_nifti);

%% PREP ROI MASK
% use custom function that loads brainstemnavigator ROI for Bstem_SC
mask = select_this_atlas_subset(region);

%% LOAD FMRI TIMESERIES FOR A SINGLE SUBJECT/TASK, but ACROSS RUN

% the preprocessing must be done on the whole-brain data because, e.g., the whole-image RMSSD flagger assumes whole brain data
preprocessed_dat = preproc_fmri_data_canlabtools(trs_to_use, n_runs, tr_duration, paths_nifti, paths_confounds);
masked_dat = apply_mask(preprocessed_dat, mask);

%% SAVE OUT MASKED (smaller) DATA
% these are basically only ever getting read in by other canlabtools analysis scripts I think
% so fine to write them as .mat
% pull out just the timeseries from the fmri_data object
% and transpose it to get time on the row and voxel on the column
DATA = masked_dat.dat';
save(out_path, 'DATA')
