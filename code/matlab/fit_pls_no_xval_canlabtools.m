%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO ESTIMATE THE BETAS ACROSS THE WHOLE DATASET FOR A SINGLE ENCODING MODEL
% ONLY FOR USE TO LOOK AT SPATIAL PATTERNS OF BETAS. NOT TO EVALUATE MODEL PERFORMANCE!!!
% expects inputs with multiple runs per subject
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
% out_path_betas: for the model beta matrix (time x voxels).
% recommended to put the ROI in the file name as well so that when this is run
% multiple times for multiple ROIs, the files will be distinct. 

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
[bold_masked_allsubjs, ~, ~] = load_fmri_data_for_pls_allsubs(paths_masked);

%% PLS THEM TOGETHER!!!
disp('Preparing to fit model')
% fit NOT CROSS-VALIDATED PLS for this ROI
% the size should be determined by this point so preallocate to be good girls
n_voxels = width(bold_masked_allsubjs);

% you can change the first value to change the number of default comps if you like
% this is set within the loop so if the current ROI has fewer voxels than the number of default comps
% then only that ROI will have fewer comps accordingly
n_pls_comps = min([100, n_voxels, width(activations)]);
fprintf('Fitting with %03d PLS components\n', n_pls_comps)

[~,~,~,~,betas] = plsregress(activations, bold_masked_allsubjs, n_pls_comps);
% aargh linear algebra: the beta matrix will have n_units (x matrix) ROWS and n_voxels (y matrix) COLUMNS

% NO HELD-OUT DATA FOR YHATS OR PRED-OBS-CORRELATIONS BECAUSE THIS IS FOR THE BETAS ONLY!!!

%% SAVE OUT ONE SET OF BETAS (NOT BY SUBJECT REMEMBER)

writematrix(betas, out_path_betas);
