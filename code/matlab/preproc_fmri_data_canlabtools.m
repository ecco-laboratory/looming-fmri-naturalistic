function preprocessed_data = preproc_fmri_data_canlabtools(trs_to_use, n_runs, tr_duration, paths_nifti, paths_confounds)
    % this is designed to read in one subject/task's worth of fMRI data, across all of their runs
    % and then preprocess and return the whole-brain data for masking and etc elsewhere
    % NB: 3 runs of the naturalistic task takes up nearly 3 GB of memory per subject to read in
    % keep in mind when you're setting slurm memory caps to run this

    addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
    addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

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
    preprocessed_data = canlab_connectivity_preproc(bold, 'bpf', [.008 1/2], tr_duration, 'no_plots');
end
