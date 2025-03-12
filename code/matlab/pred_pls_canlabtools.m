%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% This script is designed to estimate predicted-observed correlations
% for a _combined_ encoding model pre-fit on a given ROI
% but specifically with one component model's unit predictors "lesioned"

% paths_masked: derived from target. 
% nested cell array of paths to already masked and canlabtools-preprocessed ROI timeseries. 
% one per subject, already goes across run
% tr_duration: in seconds, for spm_hrf()
% paths_activations_1 OR paths_activations_2: derived from target. 
% nested cell array to csvs of 2D matrices already interpolated into TR frequency
% nested FIRST by subject, THEN by run.
% whichever one is set will be the one that is NOT lesioned. 
% whichever one is left unset will be created as all 0s and thus lesioned for prediction.
% path_betas: derived from target. 
% single path to table of cross-validated PLS betas fit to this ROI for a given encoding model
% from fit_pls_canlabtools. dimensions of subj-fold x encoding model unit (including intercept) on the rows,
% voxel on the column (with initial indicator column for subj-fold)
% the source encoding model should be a joint model containing the activations specified in paths_activations_1 or 2
% 1 vs 2 matters: must set it in the location corresponding to the betas of the joint model
% out_path_perf: for the performance on held out subjects.
% out_path_pred: for the model-predicted multivoxel timecourses.
% recommended to put the ROI in the file names as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 

n_subjs = length(paths_masked);

%% LOAD PREVIOUSLY FIT PLS BETAS
% this script does not check that the betas were fit on the correct ROI
% that should be managed in the targets wrapper
betas = readmatrix(path_betas);

% pull the indicator column off into its own object
beta_subj_indices = betas(:,1);
betas = betas(:, 2:end);
% when calculating the width of the original encoding model, don't count the intercept beta
n_pls_units = height(betas(beta_subj_indices==1, :)) - 1;

%% LOAD/CONCATENATE/HRF-CONVOLVE ENCODING MODEL TIMECOURSES
% yeah you do need to do this again here because you need the activations read in to multiply against the pre-fit betas lol
% ACTIVATIONS HAVE ALREADY BEEN INTERPOLATED INTO TR FREQUENCY 
% BY THE TIME THEY GET HERE, BUT NOT CONVOLVED YET
% so keep the dims separate by run until you convolve

if exist('paths_activations_1', 'var') == 1
    activations = load_encoding_activations_allsubs(paths_activations_1, tr_duration);
    % paste on placeholder 0 (lesioned) activations fitting the dimensions of missing activations_2 to the RIGHT
    activations = [activations, zeros(height(activations), n_pls_units-width(activations))];
elseif exist('paths_activations_2', 'var') == 1
    activations = load_encoding_activations_allsubs(paths_activations_2, tr_duration);
    % in this case, paste on placeholder activations fitting the dimensions of missing activations_1 to the LEFT
    activations = [zeros(height(activations), n_pls_units-width(activations)), activations];
end

%% LOAD PRE-PREPROCESSED AND MASKED FMRI TIMESERIES FOR THIS ROI
% SO IT IS HIGHLY INCUMBENT ON YOU TO MAKE SURE THAT THE ACTIVATIONS AND BOLDS COME IN IN 
% MATCHED ORDER BECAUSE PAST THIS POINT YOU ARE JUST ASSUMING THESE INDICES LINE UP
disp('Loading pre-masked and pre-preprocessed timeseries')
% fencepost for first subject to allow preallocation for all subjects
bold = load(paths_masked{1});
% assume this is going to be true for every subject. same timeseries duration
bold_height = height(bold.DATA);

bold_masked_allsubjs = zeros(n_subjs*bold_height, width(bold.DATA));
subj_indices = zeros(n_subjs*bold_height, 1);
this_subj_indices = 1:bold_height;
bold_masked_allsubjs(this_subj_indices, :) = bold.DATA;
subj_indices(this_subj_indices, 1) = 1;

% loop over the rest of the subjects
for i=2:n_subjs
    clear bold
    bold = load(paths_masked{i});
    % whatever. I doubt failure
    this_subj_indices = (1:bold_height) + (bold_height * (i-1));
    bold_masked_allsubjs(this_subj_indices, :) = bold.DATA;
    subj_indices(this_subj_indices, 1) = i;
    
end
clear bold this_subj_indices

%% GET PLS PREDS FOR EACH HELD-OUT SUBJECT
% for each encoding model, fit the PLS

disp('Preparing to calculate model predictions')
% the size should be determined by this point so preallocate to be good girls
n_voxels = width(bold_masked_allsubjs);

yhat = zeros(size(bold_masked_allsubjs));
pred_obs_corr = zeros(n_subjs, n_voxels);

fprintf('Current held-out subject:           ')
for k=1:n_subjs
    fprintf('\b\b\b\b\b\b\b\b\b\b%03d of %03d', k, n_subjs)
    % for everything with timecourse data. the betas have their own subj indices
    test_idx = subj_indices==k;

    % only need to calc yhat here because the betas have already been fit :3
    % again! one of the input encoding models' columns have been zeroed out
    yhat(test_idx,:) = [ones(height(activations(test_idx,:)), 1), activations(test_idx,:)] * betas(beta_subj_indices==k, :);
    
    % has to be like this bc corr(X, Y) correlates each pair of columns (here, voxels)
    % but we only care about correlating each voxel's real data to its own predicted data
    % so diag() keeps only each voxel to itself
    pred_obs_corr(k,:) = diag(corr(yhat(test_idx,:), bold_masked_allsubjs(test_idx,:)));
    
end
fprintf('\n')

%% WRITE OUT

writematrix(pred_obs_corr, out_path_perf);
writematrix([subj_indices, yhat], out_path_pred);
