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
% region: in canlabcore atlas syntax
% out_path: for the performance on held out subjects

n_subjs = length(paths_nifti);

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

% flag timepoints to exclude in a binary variable
exclude = true(max(trs_to_use),1); 
exclude(trs_to_use)=false; 
exclude = repmat(exclude,3,1);

for i=1:n_subjs
    
    % append the volume indices for the non-discarded volumes here
    % using the same logic as in specify_estimate_level1
    % paths_nifti_formatted = {};
    % for j=1:length(paths_nifti{i})
    %     paths_nifti_formatted = [
    %         paths_nifti_formatted;
    %         cellstr(append(paths_nifti{i}{j}, ',', string(trs_to_use)))'
    %     ];
    %     % create the constant idx of subj number here so that it lines up with TR counts
    %     % be mindful that these are subj idx, not the actual subj id for the filenames
    %     subj_indices = [subj_indices; repmat(i, length(trs_to_use), 1)];
    % end



    % now read in just this subject's whole brain data, but stacked across runs
    % NB: 3 runs of the naturalistic task takes up nearly 3 GB of memory per subject to read in
    % which is why this is looped, to keep the total memory ceiling down 
    % by keeping/concatenating only the data from the (smaller) ROI of interest
    % NB2: fmri_data DOES accept 4D paths in SPM syntax, i.e. 'timeseries.nii,1'
    % it will throw a 'Cannot find file:' warning but DO NOT FEAR!!!
    % bold = fmri_data(paths_nifti_formatted);
    
    bold = fmri_data([paths_nifti{i}(:)]);
    bold.dat(:,exclude) = [];

    % 2024-12-18 Leaving the call here in case you want to filter later BUT don't do it right now
    % Especially don't band pass using the flynet1 studyforrest retinotopy parameters lol.
    % don't filter out the task signal!!!
    % bold = canlab_connectivity_preproc(bold, 'bpf', [.667/32 2/32],2);
    
    % masky mask
    bold_masked = apply_mask(bold, mask);
    % only now, after masking, do we append to the everybody data
    % as usual, transpose to get time on the rows and voxel on the columns
    % currently not preallocating this (yeah I know sorry) bc looping over subjects isn't that bad
    bold_masked_allsubjs = [bold_masked_allsubjs; bold_masked.dat'];
end
clear bold bold_masked

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
