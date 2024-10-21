%-----------------------------------------------------------------------
% Job saved on 09-Feb-2024 16:29:35 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%% load libraries, I think
addpath('/home/data/eccolab/Code/GitHub/spm12');
% must wake up spm shit manually bc calling through cmd line
spm('defaults','fmri');
spm_jobman('initcfg');

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% model_path

%% SPM BATCH MUMBO JUMBO
% '/home/data/eccolab/SPLaT/models/task-naturalistic/acq-mb8/sub-9901/model-endspike/smoothed-2mm/SPM.mat'
matlabbatch{1}.spm.util.voi.spmmat = {model_path};
% todo: set an f-contrast for all the task conditions so this will regress out motion etc
matlabbatch{1}.spm.util.voi.adjust = 1;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'sc';
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-4.76 -32.66 -5.1];
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 4;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [4.76 -32.66 -5.1];
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 4;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1 | i2';

%% actually execute
spm_jobman('run',matlabbatch);
