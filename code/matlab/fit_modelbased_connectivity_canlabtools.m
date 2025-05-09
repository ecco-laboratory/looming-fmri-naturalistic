%% setup
% library things
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));


%% variables that must be set in targets before the script is called
% THIS SCRIPT IS DESIGNED TO CALCULATE MODEL-BASED FUNCTIONAL CONNECTIVITY FOR A **SINGLE** HELD-OUT SUBJECT
% It takes in cross-validated model-based seed timecourses that have already been fit over all subjects, 
% so this script only needs to read in one subject at a time and train-test split separation can still be maintained. 
% which is GOOD because whew mama, the whole brain fmri_data object gets untenably large when reading in more than a few subjects at once.
% a downstream target will combine subjects' connectivity maps for t-value estimation.
% subj_num: integer of subject ID. **NOT** the cross-validation fold number (so that it will be robust to potential changes in fold structure)
% trs_to_use: vector of valid TR indices for SPM indexing into the 4D volume
% tr_duration: in seconds, for the band-pass filter preproc
% paths_nifti: derived from target. 
% cell array of paths to one subject's smoothed timeseries niftis
% paths_confounds: derived from target. 
% cell array to txts of SPM-compatible 2D matrices of confound regressors of interest
% already had the initial volumes of each removed
% BOTH paths_nifti and paths_confounds here must contain SINGLE SUBJECT DATA: nested by RUN
% seed_yhat_path: path to tabular data file containing cross-validated model-predicted multivoxel timecourses from seed region
% this will be for ALL SUBJECTS and then this script will filter it for the applicable subject
% out_path: full csv path to save the subj x whole-brain correlations out

n_runs = length(paths_nifti);

%% LOAD IN DATA FROM TARGETS PATHS
% model-predicted multivoxel timecourse for a given seed region and encoding model
% THE SEED ROI AND ENCODING MODEL ARE IMPLICITLY BY THE PATH TO THE TIMECOURSE
% this is the output written to out_path_pred by fit_pls_canlabtools, written out for _all subjects_
% so you will need to read in the same file every time for each connectivity subject but that's fine bc it's fast

% 2025-02-18: the first column of this is now subject indices
seed_yhat = readmatrix(seed_yhat_path);
% average the predicted response across "voxels" of the seed region 
% to get a sub-time vector for the whole seed region
seed_yhat_mean = squeeze(mean(seed_yhat(:, 2:end),2));

%% READ AND CANLABTOOLS-PREPROCESS EVERYONE'S WHOLE BRAIN DATA :(

% fencepost for first subject so we can preallocate the whole-brain matrix
bold_preprocessed = preproc_fmri_data_canlabtools(trs_to_use, n_runs, tr_duration, paths_nifti, paths_confounds);
dat_preprocessed = bold_preprocessed.dat';
dat_preprocessed = zscore(dat_preprocessed);

subj_indices = seed_yhat(:, 1); % may as well since it already exists

%% CORRELATE THIS SUBJECT'S MODEL-SEED TIMECOURSE WITH THEIR WHOLE-BRAIN DATA
fprintf('Correlating whole-brain data with model-seed timecourse...\n')
% first pull the subj number from the BOLD path
[~, file_nifti, ~] = fileparts(paths_nifti{1});
% the substring positions assume it starts with the prefix "smoothed_4mm_"
this_subj_num = str2double(extractBetween(convertCharsToStrings(file_nifti), 18, 21));

% need to index seed_yhat_mean correctly to get just the relevant subject's model-predicted seed timecourse
% based on matrix dimensions, this correlation will come out as a row vector with voxels on the columns
vox_corr_wb = corr(seed_yhat_mean(subj_indices==this_subj_num), dat_preprocessed);
% NB! transpose to column vector for writing out because R packages will read this in a little better later
vox_corr_wb = vox_corr_wb';

%% write out whole-brain correlations. to tabular data though
% do the between model diffs in another script. this operates on one encoding model at a time
fprintf('Saving results to file!\n')
writematrix(vox_corr_wb, out_path);
