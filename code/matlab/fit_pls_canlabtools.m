%% READ ME! LOGIC EXPLANATION
% This script is written in a funky-ish way, because the nifti-reading chunk takes a while (~2 min/subject).
% It doesn't seem like a huge time save to write the fmri_data objects to .mat 
% and read them in that way, because they're still huge in memory.
% Because of that, we want to read the niftis in as _infrequently_ as possible.
% Depending on the situation, we may have some arbitrary number of encoding model activations
% and arbitrary number of brain ROIs within which we want to model 
% the BOLD timeseries by each encoding model one-by-one.
% This script accepts an arbitrary set of activation timeseries from one or more encoding models,
% and an arbitrary number of ROIs,
% reads in the BOLD data _once,_ and then fits each encoding model to each ROI's data one by one.
% Ordinarily such a script would read in one encoding model, read in BOLD, 
% mask for one ROI, and then fit, with repetition in the wrapper call. 
% but for the reasons I described, 
% we need to not re-read the BOLD data in for each encoding model/ROI.

%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS MULTIPLE SUBJECTS AND ENCODING MODELS!!!
% with multiple runs per subject
% the nested cell arrays are first by subject, then by run within subject
% each of these should come in from targets calling script

% paths_masked: derived from target. 
% nested cell array of paths to already masked and canlabtools-preprocessed ROI timeseries. 
% one per subject, already goes across run
% tr_duration: in seconds, for spm_hrf()
% paths_activations: derived from target. 
% nested cell array to csvs of 2D matrices already interpolated into TR frequency
% THIS one is nested slightly differently: FIRST by encoding model, THEN by subject, THEN by run.
% out_paths: for the performance on held out subjects. cell array by encoding model.
% will write one file per encoding model. recommended to put the ROI in the file names as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 
% if this variable does not exist, performance will not be written out to file.

n_subjs = length(paths_masked);
n_runs = length(paths_activations{1});

%% LOAD/CONCATENATE/HRF-CONVOLVE ENCODING MODEL TIMECOURSES
% ACTIVATIONS HAVE ALREADY BEEN INTERPOLATED INTO TR FREQUENCY 
% BY THE TIME THEY GET HERE, BUT NOT CONVOLVED YET
% so keep the dims separate by run until you convolve

activations = [];
disp('Loading encoding model activation timecourses')
% then subject
for j=1:length(paths_activations)
    fprintf('Current subject: %03d of %03d\n', j, n_subjs)
    % then run
    for k=1:length(paths_activations{j})
        activations_this_run = readmatrix(paths_activations{j}{k});
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

% TODO: decide whether this one is single subject or looping over subjects
% MT is leaning toward doing this one single subject and having a loop called elsewhere

%% LOOP OVER ROIS AND PLS THEM TOGETHER!!!
% for each encoding model, fit the PLS

disp('Preparing to fit model')
% fit cross-validated PLS for this ROI
% the size should be determined by this point so preallocate to be good girls
n_voxels = width(bold_masked_allsubjs);
yhat = zeros(size(bold_masked_allsubjs));
pred_obs_corr = zeros(n_subjs, n_voxels);

% you can change the first value to change the number of default comps if you like
% this is set within the loop so if the current ROI has fewer voxels than the number of default comps
% then only that ROI will have fewer comps accordingly
n_pls_comps = min(100, int8(n_voxels));
fprintf('Fitting with %03d PLS components\n', n_pls_comps)

for k=1:n_subjs
    fprintf('Current held-out subject: %03d of %03d\n', k, n_subjs)
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
if exist('out_path', 'var') == 1
    writematrix(pred_obs_corr, out_path);
end

clear mask
