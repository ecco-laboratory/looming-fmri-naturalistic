%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS MULTIPLE SUBJECTS!!!
% with multiple runs per subject
% the nested cell arrays are first by subject, then by run within subject
% each of these should come in from targets calling script

% paths_nifti: derived from target. 
% nested cell array of paths to smoothed timeseries niftis
% trs_to_use: vector of valid TR indices for SPM indexing into the 4D volume
% tr_duration: in seconds, for spm_hrf()
% paths_activations: derived from target. 
% nested cell array to csvs of 2D matrices already interpolated into TR frequency
% paths_confounds: derived from target. 
% nested cell array to txts of SPM-compatible 2D matrices of confound regressors of interest
% already had the initial volumes of each removed
% region: in canlabcore atlas syntax
% out_path: for the performance on held out subjects

n_subjs = length(paths_nifti);
n_runs = length(paths_nifti{1});

%% LOAD/CONCATENATE/HRF-CONVOLVE ENCODING MODEL TIMECOURSES

% ALREADY BEEN INTERPOLATED INTO TR FREQUENCY BY THE TIME THEY GET HERE
% BUT NOT CONVOLVED YET
% so keep the dims separate by run until you convolve
activations = [];
for i=1:length(paths_activations)
    for j=1:length(paths_activations{i})
        activations_this_run = readmatrix(paths_activations{i}{j});
        % Convolve features with HRF here
        % holy shit... matlab anonymous functions
        conv_activations = arrayfun(@(i) conv(double(activations_this_run(:, i)), spm_hrf(tr_duration)), 1:size(activations_this_run, 2), 'UniformOutput', false);
        conv_activations = cell2mat(conv_activations);
        % trim off the tail introduced by the convolution
        conv_activations = conv_activations(1:height(activations_this_run), :);
        % now after the run activations have been convolved we can concatenate onto the main
        % should end up with a 2D array where time x subject is on the rows and unit is on the cols
        activations = [activations; conv_activations];
    end
end
clear activations_this_run conv_activations

%% LOAD FMRI TIMESERIES

% SO IT IS HIGHLY INCUMBENT ON YOU TO MAKE SURE THAT THE ACTIVATIONS AND BOLDS COME IN IN 
% MATCHED ORDER BECAUSE PAST THIS POINT YOU ARE JUST ASSUMING THESE INDICES LINE UP
subj_indices = [];
bold_masked_allsubjs = [];
mask = select_atlas_subset(load_atlas('canlab2018'), {region});

% flag timepoints to exclude by logical indexing
exclude = true(max(trs_to_use),1); 
exclude(trs_to_use) = false; 
% this should handle discarding the first several volumes for each run independently
exclude = repmat(exclude, n_runs, 1);

for i=1:n_subjs
    % now read in just this subject's whole brain data, but stacked across runs
    % NB: 3 runs of the naturalistic task takes up nearly 3 GB of memory per subject to read in
    % which is why this is looped, to keep the total memory ceiling down 
    % by keeping/concatenating only the data from the (smaller) ROI of interest
    
    bold = fmri_data([paths_nifti{i}(:)]);
    bold.dat(:,exclude) = [];
    % create the constant idx of subj number here so that it lines up with TR counts
    % be mindful that these are subj idx, not the actual subj id for the filenames
    subj_indices = [subj_indices; repmat(i, width(bold.dat), 1)];
    
    % masky mask
    bold_masked = apply_mask(bold, mask);

    % light preprocessing

    % first, use canlabtools method to ID volumes with a big sequential jump in RMSSD
    % the method generates a movie for interactive viewing by default. hence turning it off in the args
    [rmssd, rmssd_outlier_regressor_matrix] = rmssd_movie(bold_masked,'showmovie',false,'nodisplay');
    
    % then read in and row-run-stack fmriprep motion regressors
    confounds = []; session_means =[];
    for j=1:n_runs
        confounds = [confounds; readmatrix(paths_confounds{i}{j})];
        session_means = [session_means; j*ones(height(readmatrix(paths_confounds{i}{j})),1)];
    end

    % append RMSSD regressors to fmriprep motion regressors, attach to the fmri_data obj, regress
    % 2024-12-18: note that we aren't doing any X-pass filtering
    bold_masked.covariates =[confounds, rmssd_outlier_regressor_matrix, condf2indic(session_means)];
    [preprocessed_dat] = canlab_connectivity_preproc(bold_masked,'bpf', [.008 1/8],tr_duration, 'no_plots');

    % only now, after masking and preprocessing, do we append to the everybody data
    % as usual, transpose to get time on the rows and voxel on the columns
    % currently not preallocating this (yeah I know sorry) bc looping over subjects isn't that bad
    bold_masked_allsubjs = [bold_masked_allsubjs; preprocessed_dat.dat'];
end
clear bold bold_masked preprocessed_dat

% TODO: decide whether this one is single subject or looping over subjects
% MT is leaning toward doing this one single subject and having a loop called elsewhere

%% PLS THEM TOGETHER!!!
% you can change the first value to change the number of default comps if you like
n_pls_comps = min(20, size(bold_masked_allsubjs,2));

% the size should be determined by this point so preallocate to be good girls
n_voxels = width(bold_masked_allsubjs);
yhat = zeros(size(bold_masked_allsubjs));
pred_obs_corr = zeros(n_subjs, n_voxels);

for k=1:n_subjs
    train_idx = subj_indices~=k;
    test_idx = ~train_idx;
    [~,~,~,~,beta_cv] = plsregress(activations(train_idx,:), bold_masked_allsubjs(train_idx,:), n_pls_comps);
    
    % ones() prepends an intercept column to the selected chunk of conv_features to be used for predicting against the betas
    yhat(test_idx,:)=[ones(sum(test_idx),1) activations(test_idx,:)]*beta_cv;
    
    % has to be like this bc corr(X, Y) correlates each pair of columns (here, voxels)
    % but we only care about correlating each voxel's real data to its own predicted data
    % so diag() keeps only each voxel to itself
    pred_obs_corr(k,:)=diag(corr(yhat(test_idx,:), bold_masked_allsubjs(test_idx,:)));
    
end

%% SAVE OUT RELEVANT RESULTS
writematrix(pred_obs_corr, out_path);
