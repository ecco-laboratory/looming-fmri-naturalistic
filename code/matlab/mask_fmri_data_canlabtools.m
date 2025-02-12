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
% NB: 3 runs of the naturalistic task takes up nearly 3 GB of memory per subject to read in
% keep in mind when you're setting slurm memory caps to run this

% flag discard acquisition timepoints to exclude by logical indexing
exclude = true(max(trs_to_use),1); 
exclude(trs_to_use) = false; 
% this should handle discarding the first several volumes for each run independently
exclude = repmat(exclude, n_runs, 1);
% load the actual niftis (this takes a little bit!)
bold = fmri_data([paths_nifti(:)]);
% discard those timepoints
bold.dat(:,exclude) = [];

%% PREPROCESS WITH CANLABTOOLS
% first, use canlabtools method to ID volumes with a big sequential jump in RMSSD
% the method generates a movie for interactive viewing by default. hence turning it off in the args
[~, rmssd_outlier_regressor_matrix] = rmssd_movie(bold,'showmovie',false,'nodisplay');

% then read in and row-run-stack fmriprep motion regressors
confounds = []; session_means = [];
for i=1:n_runs
    these_confounds = readmatrix(paths_confounds{i});
    confounds = [confounds; these_confounds];
    session_means = [session_means; repmat(i, height(these_confounds), 1)];
end

% append RMSSD regressors to fmriprep motion regressors, attach to the fmri_data obj, regress
bold.covariates =[confounds, rmssd_outlier_regressor_matrix, condf2indic(session_means)];
% 2025-01: also doing some band-pass filtering around the range of expected task activation
preprocessed_dat = canlab_connectivity_preproc(bold, 'bpf', [.008 1/2], tr_duration, 'no_plots');

% masky mask
masked_dat = apply_mask(preprocessed_dat, mask);

%% SAVE OUT MASKED (smaller) DATA
% these are basically only ever getting read in by other canlabtools analysis scripts I think
% so fine to write them as .mat
% pull out just the timeseries from the fmri_data object
% and transpose it to get time on the row and voxel on the column
DATA = masked_dat.dat';
save(out_path, 'DATA')
