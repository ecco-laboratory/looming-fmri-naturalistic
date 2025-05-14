%-----------------------------------------------------------------------
% Job saved on 14-Nov-2023 15:39:17 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%% load libraries, I think
addpath('/home/data/eccolab/Code/GitHub/spm12');
% must wake up spm shit manually bc calling through cmd line
spm('defaults','fmri');
spm_jobman('initcfg');

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% out_path for the SPM.mat (just the folder)
% paths_nifti (is derived from a target)
% trs_to_use, already a vector of valid tr nums
% paths_onsets (is derived from a target)
% paths_confounds (is derived from a target), already subset by columns and with first 8s-ish TRs' worth of confounds removed
% basis

%% SPM BATCH #1 MUMBO JUMBO
% specify ze main glm
% SPM.mat gets written into this dir (and you can't change the filename) so make sure the dir is specifically pointing to one sequence's folder
batch_spec_est{1}.spm.stats.fmri_spec.dir = {out_path};
batch_spec_est{1}.spm.stats.fmri_spec.timing.units = 'secs';
batch_spec_est{1}.spm.stats.fmri_spec.timing.RT = tr_duration; % this is TR but SPM naming is flipped like a silly goose
batch_spec_est{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
batch_spec_est{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% loop over runs, popping in the path for each run's niftis/onsets/confound regressors
for i=1:length(paths_confounds)
    % format with the tr_nums here
    % so that the input call string from R doesn't get CRAZY long
    paths_nifti_formatted = cellstr(append(paths_nifti{i}, ',', string(trs_to_use)));
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).scans = paths_nifti_formatted;
    %%
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    % set by targets in the calling script
    % batch_spec_est{1}.spm.stats.fmri_spec.sess.multi = cellstr(path_onsets);
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).multi = paths_onsets(i);
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    % set by targets in the calling script
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).multi_reg = paths_confounds(i);
    batch_spec_est{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
end

batch_spec_est{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
if strcmp(basis, 'fir')
    batch_spec_est{1}.spm.stats.fmri_spec.bases.fir.length = 16;
    batch_spec_est{1}.spm.stats.fmri_spec.bases.fir.order = 16;
else
    batch_spec_est{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end
batch_spec_est{1}.spm.stats.fmri_spec.volt = 1;
batch_spec_est{1}.spm.stats.fmri_spec.global = 'None';
batch_spec_est{1}.spm.stats.fmri_spec.mthresh = 0.8;
batch_spec_est{1}.spm.stats.fmri_spec.mask = {''};
batch_spec_est{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% estimate ze model
batch_spec_est{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
batch_spec_est{2}.spm.stats.fmri_est.write_residuals = 0;
batch_spec_est{2}.spm.stats.fmri_est.method.Classical = 1;

%% execute the first batch
spm_jobman('run',batch_spec_est);

%% CONSTRUCT THE CONTRASTS FROM THE DESIGN MATRIX
% we load the design matrix in SPM.mat into matlab (not within SPM)
% so that we can use the condition names already set in the matrix to set contrasts
% this is already downstream of combining across runs so that's handled!
% the first batch is defined and run first so that the SPM.mat exists by this point
load(fullfile(out_path, 'SPM.mat'));
designmat_names = convertCharsToStrings(SPM.xX.name);

% for the controlled task
% important that these are strings defined with double quotes
% currently not doing any f-contrasts 
% but if you were it would have to be a squarish matrix with a 1 for each stim condition

% Also, make contrasts for loom x animal vs baseline to average over spatial location for testing naturalistic encoding models on controlled betas
if contains(out_path, "controlled")
    % important that these are strings defined with double quotes
    consess = calc_contrasts_level1_controlled(designmat_names);
elseif contains(out_path, "naturalistic")
    consess = calc_contrasts_level1_naturalistic(designmat_names);
end
% set the sessrep field as 'none' for all of them
% because the same exact conditions don't appear in every run
% so we can't use the replicate across session setting
for contrast_num = 1:length(consess)
    consess{contrast_num}.tcon.sessrep = 'none';
end

%% SPM BATCH #2 MUMBO JUMBO
% specify ze contrasts
% when it was all in the same matlabbatch and calling from dependency: cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
% we can't do that anymore because we need this batch to run after the other one is finished
% so that we can read in the design matrix ourselves
batch_contrast{1}.spm.stats.con.spmmat(1) = cellstr(fullfile(out_path, 'SPM.mat'));
batch_contrast{1}.spm.stats.con.consess = consess;
batch_contrast{1}.spm.stats.con.delete = 1;

%% execute the second batch
spm_jobman('run',batch_contrast);
