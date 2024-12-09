%-----------------------------------------------------------------------
% Job saved on 31-Jan-2024 16:03:11 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%% load libraries, I think
addpath('/home/data/eccolab/Code/GitHub/spm12');
% must wake up spm shit manually bc calling through cmd line
spm('defaults','fmri');
spm_jobman('initcfg');
%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% paths_nifti: make sure the it comes in as a cell array with the `,n` TR digits appended to each
% kernel: width in mm. this script will set an isotropic kernel
% out_prefix: set in the targets file so that it's forcibly consistent
%% SPM BATCH MUMBO JUMBO
matlabbatch{1}.spm.spatial.smooth.data = paths_nifti;
matlabbatch{1}.spm.spatial.smooth.fwhm = repelem(kernel, 3); % in mm, x y z
matlabbatch{1}.spm.spatial.smooth.dtype = 0; % return with unchanged data type 
matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit masking
matlabbatch{1}.spm.spatial.smooth.prefix = out_prefix; % why god couldn't they have done a suffix

%% actually execute
spm_jobman('run',matlabbatch);
