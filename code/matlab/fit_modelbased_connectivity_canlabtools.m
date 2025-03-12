%% setup
% library things
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

% this script ONLY processes studyforrest data
% and the SPM preprocessing script upstream of this will only output the realigned/normalized niftis
% back to the dir where the original niftis are saved
% so we are left going into the studyforrest folder
studyforrest_dir = '../ignore/datasets/studyforrest-data-phase2';

%% variables that must be set in targets before the script is called
% trs_to_use: vector of valid TR indices for SPM indexing into the 4D volume
% tr_duration: in seconds, for the band-pass filter preproc
% paths_nifti: derived from target. 
% nested cell array of paths to smoothed timeseries niftis
% paths_confounds: derived from target. 
% nested cell array to txts of SPM-compatible 2D matrices of confound regressors of interest
% already had the initial volumes of each removed
% BOTH paths_nifti and paths_confounds here must contain ALL SUBJECTS DATA: nested FIRST by SUBJECT and THEN by RUN
% seed_yhat_path: path to tabular data file containing cross-validated model-predicted multivoxel timecourses from seed region
% out_path: full csv path to save the subj x whole-brain correlations out

n_subjs = length(paths_nifti);
n_runs = length(paths_nifti{1});

%% EMERGENCY SETTING THOSE VARIABLES IF NOT RUNNING FROM TARGETS
% sc_data_path = fullfile(studyforrest_dir, 'fmri_data_canlabtooled_sc.mat');
% studyforrest_activation_path = fullfile(studyforrest_dir, 'flynet_convolved_timecourses.csv');
% out_fstring = '%smap_flynet_connectivity_contrast.nii';

%% LOAD IN DATA FROM TARGETS PATHS
% model-predicted multivoxel timecourse for a given seed region and encoding model
% THE SEED ROI AND ENCODING MODEL ARE IMPLICITLY BY THE PATH TO THE TIMECOURSE
% 2025-02-18: the first column of this is now subject indices
seed_yhat = readmatrix(seed_yhat_path);
% average the predicted response across "voxels" of the seed region 
% to get a sub-time vector for the whole seed region
seed_yhat_mean = squeeze(mean(seed_yhat(:, 2:end),2));

%% READ AND CANLABTOOLS-PREPROCESS EVERYONE'S WHOLE BRAIN DATA :(

% fencepost for first subject so we can preallocate the whole-brain matrix
bold_preprocessed = preproc_fmri_data_canlabtools(trs_to_use, n_runs, tr_duration, paths_nifti{1}, paths_confounds{1});
dat_preprocessed = bold_preprocessed.dat';
bold_height = height(dat_preprocessed);

WB_DATA = zeros(n_subjs*bold_height, width(dat_preprocessed));

subj_indices = seed_yhat(:, 1); % may as well since it already exists
WB_DATA(subj_indices==1, :) = zscore(dat_preprocessed);

% loop over the rest of the subjects
for s=2:n_subjs
    clear bold_preprocessed dat_preprocessed
    bold_preprocessed = preproc_fmri_data_canlabtools(trs_to_use, n_runs, tr_duration, paths_nifti{s}, paths_confounds{s});
    dat_preprocessed = bold_preprocessed.dat';
    % no masking. whole brain. not parcel-averaged. all voxels. huge.
    WB_DATA(subj_indices==s, :) = zscore(dat_preprocessed);
end
clear bold_preprocessed dat_preprocessed

%% CORRELATE MODEL-SEED TIMECOURSE WITH WHOLE-BRAIN DATA, BY SUBJECT

% leave-one-subject-out cross-validated functional connectivity. great
vox_corr_wb = zeros(n_subjs, width(WB_DATA));

fprintf('Correlating whole-brain data with model-seed timecourse...\n')
for s = 1:n_subjs
    % correlate the model-predicted seed timecourse with each whole-brain voxel timecourse
    % you have to loop over subjects to get the prediction for each subject. 
    % if you just corr on the whole thing it gives you the seed vox-target vox correlation across all subjects
    % we now need to save the mean y-hats out for the permutation testing later as well
    vox_corr_wb(s,:) = corr(seed_yhat_mean(subj_indices==s), WB_DATA(subj_indices==s,:));
end

%% write out whole-brain correlations. to tabular data though
% do the between model diffs in another script. this operates on one encoding model at a time
fprintf('Saving results to file!\n')
writematrix(vox_corr_wb, out_path);
