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
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

% this atlas is not implemented through CanlabCore/Neuroimaging_Pattern_Masks 
% so we're just gonna have to literally read in some niftis from this folder
BrainstemNavigator_path = '/home/data/shared/BrainstemNavigator/0.9/2a.BrainstemNucleiAtlas_MNI/labels_thresholded_binary_0.35';

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS MULTIPLE SUBJECTS AND ENCODING MODELS!!!
% with multiple runs per subject
% the nested cell arrays are first by subject, then by run within subject
% each of these should come in from targets calling script

% paths_nifti: derived from target. 
% nested cell array of paths to smoothed timeseries niftis
% trs_to_use: vector of valid TR indices for SPM indexing into the 4D volume
% tr_duration: in seconds, for spm_hrf()
% paths_activations: derived from target. 
% nested cell array to csvs of 2D matrices already interpolated into TR frequency
% THIS one is nested slightly differently: FIRST by encoding model, THEN by subject, THEN by run.
% paths_confounds: derived from target. 
% nested cell array to txts of SPM-compatible 2D matrices of confound regressors of interest
% already had the initial volumes of each removed
% regions: cell array in canlabcore atlas syntax. Each cell will be treated as a separate mask
% and the analysis will be repeated for each mask. 
% this is to reduce the number of times the full dataset needs to be read in from nifti
% out_paths: for the performance on held out subjects. IMPORTANT: nested cell array, first by encoding model and then by ROI.
% will write one file per encoding model x ROI. 
% if this variable does not exist, performance will not be written out to file.

n_subjs = length(paths_nifti);
n_runs = length(paths_nifti{1});

%% LOAD/CONCATENATE/HRF-CONVOLVE ENCODING MODEL TIMECOURSES
% ACTIVATIONS HAVE ALREADY BEEN INTERPOLATED INTO TR FREQUENCY 
% BY THE TIME THEY GET HERE, BUT NOT CONVOLVED YET
% so keep the dims separate by run until you convolve

activations = cell(size(paths_activations));
% first loop over encoding model
for i=1:length(paths_activations)
    fprintf('Loading activations for model: %02d of %02d\n', i, length(paths_activations))
    % then subject
    for j=1:length(paths_activations{i})
        fprintf('Current subject: %03d of %03d\n', j, n_subjs)
        % then run
        for k=1:length(paths_activations{i}{j})
            activations_this_run = readmatrix(paths_activations{i}{j}{k});
            % Convolve features with HRF here
            % holy shit... matlab anonymous functions
            conv_activations = arrayfun(@(i) conv(double(activations_this_run(:, i)), spm_hrf(tr_duration)), 1:size(activations_this_run, 2), 'UniformOutput', false);
            conv_activations = cell2mat(conv_activations);
            % trim off the tail introduced by the convolution
            conv_activations = conv_activations(1:height(activations_this_run), :);
            % now after the run activations have been convolved we can concatenate onto the main
            % should end up with a 2D array where time x subject is on the rows and unit is on the cols
            activations{i} = [activations{i}; conv_activations];
        end
    end
end
clear activations_this_run conv_activations

%% PREP ROI MASKS

masks = cell(length(regions), 1);
for i=1:length(regions)
    % First, load the mask for this ROI
    % for superior colliculus, manually use BrainstemNavigator ROI instead
    if any(strcmp(regions{i},'Bstem_SC'))
        mask = fmri_data(fullfile(BrainstemNavigator_path, 'SC_l.nii'));
        mask_r = fmri_data(fullfile(BrainstemNavigator_path, 'SC_r.nii'));
        mask.dat = mask.dat + mask_r.dat;
        % and then explicitly exclude PAG from that SC ROI
        pag = load_atlas('Kragel2019PAG');
        pag = resample_space(pag,mask);
        mask.dat(pag.dat>0) = 0;
        clear mask_r pag
    else
        mask = select_atlas_subset(load_atlas('canlab2018'), regions(i));
    end
    masks{i} = mask;
end
clear mask

%% LOAD FMRI TIMESERIES
% SO IT IS HIGHLY INCUMBENT ON YOU TO MAKE SURE THAT THE ACTIVATIONS AND BOLDS COME IN IN 
% MATCHED ORDER BECAUSE PAST THIS POINT YOU ARE JUST ASSUMING THESE INDICES LINE UP

% flag discard acquisition timepoints to exclude by logical indexing
exclude = true(max(trs_to_use),1); 
exclude(trs_to_use) = false; 
% this should handle discarding the first several volumes for each run independently
exclude = repmat(exclude, n_runs, 1);

% now read in 1 subject's whole brain data at a time, stacked across runs
% NB: 3 runs of the naturalistic task takes up nearly 3 GB of memory per subject to read in
% which is why this is looped, to keep the total memory ceiling down 
% by keeping/concatenating only each subject's data from the (smaller) ROIs of interest
subj_indices = [];
% bold_allsubjs = [];
% needs to be a cell array bc each ROI has different voxel dims
bold_masked_allsubjs = cell(length(masks), 1);

for i=1:n_subjs
    
    bold = fmri_data([paths_nifti{i}(:)]);
    bold.dat(:,exclude) = [];
    % create the constant idx of subj number here so that it lines up with TR counts
    % be mindful that these are subj idx, not the actual subj id for the filenames
    subj_indices = [subj_indices; repmat(i, width(bold.dat), 1)];
    
    % light preprocessing
    
    % first, use canlabtools method to ID volumes with a big sequential jump in RMSSD
    % the method generates a movie for interactive viewing by default. hence turning it off in the args
    [rmssd, rmssd_outlier_regressor_matrix] = rmssd_movie(bold,'showmovie',false,'nodisplay');
    
    % then read in and row-run-stack fmriprep motion regressors
    confounds = []; session_means = [];
    for j=1:n_runs
        these_confounds = readmatrix(paths_confounds{i}{j});
        confounds = [confounds; these_confounds];
        session_means = [session_means; repmat(j, height(these_confounds), 1)];
    end
    
    % append RMSSD regressors to fmriprep motion regressors, attach to the fmri_data obj, regress
    bold.covariates =[confounds, rmssd_outlier_regressor_matrix, condf2indic(session_means)];
    % 2025-01: also doing some band-pass filtering around the range of expected task activation
    preprocessed_dat = canlab_connectivity_preproc(bold, 'bpf', [.008 1/2], tr_duration, 'no_plots');
    
    % only now, after masking and preprocessing, do we append to the everybody data
    % as usual, transpose to get time on the rows and voxel on the columns
    % currently not preallocating this (yeah I know sorry) bc looping over subjects isn't that bad
    
    % WARNING! THIS IS WHOLE BRAIN DATA AND IS THUS FREAKING HUGE. COMMENTED OUT TO REDUCE MEMORY USAGE
    % bold_allsubjs = [bold_allsubjs; preprocessed_dat.dat'];
    
    % masky mask
    for j=1:length(masks)
        masked_dat = apply_mask(preprocessed_dat, masks{j});
        bold_masked_allsubjs{j} = [bold_masked_allsubjs{j}; masked_dat.dat'];
    end
    
end
clear bold bold_masked these_confounds confounds session_means rmssd rmssd_outlier_regressor_matrix preprocessed_dat masked_dat

% TODO: decide whether this one is single subject or looping over subjects
% MT is leaning toward doing this one single subject and having a loop called elsewhere

%% LOOP OVER ROIS AND PLS THEM TOGETHER!!!
% for each encoding model
for i=1:length(activations)
    % and each roi, fit the PLS
    for j=1:length(masks)
        % you can change the first value to change the number of default comps if you like
        % this is set within the loop so if the current ROI has fewer voxels than the number of default comps
        % then only that ROI will have fewer comps accordingly
        n_pls_comps = min(100, int8(size(bold_masked_allsubjs{j}, 2)));
        
        % fit cross-validated PLS for this ROI
        % the size should be determined by this point so preallocate to be good girls
        n_voxels = width(bold_masked_allsubjs{j});
        yhat = zeros(size(bold_masked_allsubjs{j}));
        pred_obs_corr = zeros(n_subjs, n_voxels);

        for k=1:n_subjs
            train_idx = subj_indices~=k;
            test_idx = ~train_idx;
            % pay attention to the counters from the outer loops: i counts encoding model, j counts ROI
            [~,~,~,~,beta_cv] = plsregress(activations{i}(train_idx,:), bold_masked_allsubjs{j}(train_idx,:), n_pls_comps);
            
            % ones() prepends an intercept column to the selected chunk of conv_features to be used for predicting against the betas
            yhat(test_idx,:)=[ones(sum(test_idx),1) activations{i}(test_idx,:)]*beta_cv;
            
            % has to be like this bc corr(X, Y) correlates each pair of columns (here, voxels)
            % but we only care about correlating each voxel's real data to its own predicted data
            % so diag() keeps only each voxel to itself
            pred_obs_corr(k,:)=diag(corr(yhat(test_idx,:), bold_masked_allsubjs{j}(test_idx,:)));
            
        end
        
        %% SAVE OUT RELEVANT RESULTS
        if exist('out_paths', 'var') == 1
            writematrix(pred_obs_corr, out_paths{i}{j});
        end
        
    end
end
clear mask
