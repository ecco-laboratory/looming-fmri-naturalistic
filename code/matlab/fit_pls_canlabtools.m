%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS MULTIPLE SUBJECTS FOR A SINGLE ENCODING MODEL
% with multiple runs per subject
% the nested cell arrays are first by subject, then by run within subject
% each of these should come in from targets calling script

% paths_masked: derived from target. 
% nested cell array of paths to already masked and canlabtools-preprocessed ROI timeseries. 
% one per subject, already goes across run
% tr_duration: in seconds, for spm_hrf()
% paths_activations_1, paths_activations_2: derived from target. 
% nested cell array to csvs of 2D matrices already interpolated into TR frequency
% nested FIRST by subject, THEN by run.
% paths_activations_1 must be provided--will fit on that one encoding model per usual.
% if paths_activations_2 also exists, will load that model's activations in as well AND THEN ROW-BIND THE PREDICTORS TOGETHER
% out_path_perf: for the performance on held out subjects.
% out_path_pred: for the model-predicted multivoxel timecourses.
% out_path_betas: for the model beta matrices.
% all will write tabular data with (held-out) subject x whatever on the rows
% recommended to put the ROI in the file names as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 
% if these variables do not exist, performance/predictions will not be written out to file.

n_subjs = length(paths_masked);

%% LOAD/CONCATENATE/HRF-CONVOLVE ENCODING MODEL TIMECOURSES
% ACTIVATIONS HAVE ALREADY BEEN INTERPOLATED INTO TR FREQUENCY 
% BY THE TIME THEY GET HERE, BUT NOT CONVOLVED YET
% so keep the dims separate by run until you convolve

activations = load_encoding_activations_allsubs(paths_activations_1, tr_duration);

if exist('paths_activations_2', 'var') == 1
    activations_2 = load_encoding_activations_allsubs(paths_activations_2, tr_duration);
    activations = [activations, activations_2];
    clear activations_2
end

%% LOAD PRE-PREPROCESSED AND MASKED FMRI TIMESERIES FOR THIS ROI
% SO IT IS HIGHLY INCUMBENT ON YOU TO MAKE SURE THAT THE ACTIVATIONS AND BOLDS COME IN IN 
% MATCHED ORDER BECAUSE PAST THIS POINT YOU ARE JUST ASSUMING THESE INDICES LINE UP
disp('Loading pre-masked and pre-preprocessed timeseries')
[bold_masked_allsubjs, subj_indices, subj_nums] = load_fmri_data_for_pls_allsubs(paths_masked);

%% PLS THEM TOGETHER!!!
disp('Preparing to fit model')
% fit cross-validated PLS for this ROI
% the size should be determined by this point so preallocate to be good girls
n_voxels = width(bold_masked_allsubjs);
n_units = width(activations);
% it will be flattened to 2d later. do this for now bc easier for me to read
% and don't forget to add 1 for the intercept :3
betas = zeros(n_units+1, n_voxels, n_subjs);
yhat = zeros(size(bold_masked_allsubjs));
pred_obs_corr = zeros(n_subjs, n_voxels);

% you can change the first value to change the number of default comps if you like
% this is set within the loop so if the current ROI has fewer voxels than the number of default comps
% then only that ROI will have fewer comps accordingly
n_pls_comps = min([100, n_voxels, width(activations)]);
fprintf('Fitting with %03d PLS components\n', n_pls_comps)

fprintf('Current held-out subject:           ')
for k=1:n_subjs
    fprintf('\b\b\b\b\b\b\b\b\b\b%03d of %03d', k, n_subjs)
    train_idx = subj_indices~=subj_nums(k);
    test_idx = ~train_idx;

    [~,~,~,~,beta_cv] = plsregress(activations(train_idx,:), bold_masked_allsubjs(train_idx,:), n_pls_comps);
    % now saving the betas out in case you want to, you know, do stuff with them
    % aargh linear algebra: the beta matrix will have n_units (x matrix) ROWS and n_voxels (y matrix) COLUMNS
    betas(:,:,k) = beta_cv;
    % ones() prepends an intercept column to the selected chunk of conv_features to be used for predicting against the betas
    % TODO: save out yhat to use as the seed timecourse for model-based connectivity analysis later
    yhat(test_idx,:) = [ones(height(activations(test_idx,:)), 1), activations(test_idx,:)]*beta_cv;
    
    % has to be like this bc corr(X, Y) correlates each pair of columns (here, voxels)
    % but we only care about correlating each voxel's real data to its own predicted data
    % so diag() keeps only each voxel to itself
    pred_obs_corr(k,:) = diag(corr(yhat(test_idx,:), bold_masked_allsubjs(test_idx,:)));
    
end
fprintf('\n')

%% SAVE OUT RELEVANT RESULTS
% 2025-02-14: It would seem that we need this to write out the yhat predictions AND the yhat-y correlations
% yhat-y correlations for "overall performance" as per usual
% but the direct yhats are needed later for, e.g., correlating predictions from different models
% and as the seed timecourses for model-based functional connectivity analyses
% the target for this will track both of these file paths together. call them independently later as you need them
if exist('out_path_perf', 'var') == 1
    writematrix(pred_obs_corr, out_path_perf);
end

% bind the subject (fold) indices to the predicted timecourses so that you can relate them to self-report etc by subject later
if exist('out_path_pred', 'var') == 1
    writematrix([subj_indices, yhat], out_path_pred);
end

% similarly, create and bind the fold indices to the cross-val betas
if exist('out_path_betas', 'var') == 1
    betas = reshape(permute(betas, [1 3 2]), [], n_voxels);
    beta_subj_indices = repelem(subj_nums, n_units+1)';
    writematrix([beta_subj_indices, betas], out_path_betas);
end
