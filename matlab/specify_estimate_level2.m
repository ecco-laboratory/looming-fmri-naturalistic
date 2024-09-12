%-----------------------------------------------------------------------
% Job saved on 22-Jul-2024 11:54:27 by cfg_util (rev $Rev: 7345 $)
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
% paths_con (is derived from a target)
% contrast_name (to put into consess below)
% the targets script handles setting the correct contrast folder name
% and specifying the con images with the correct number
% (which corresponds with the contrast order set in the 1st level spec/est batch script)

%% SPM BATCH MUMBO JUMBO
% model folder follows the format: '/home/data/eccolab/SPLaT/models/task-controlled/acq-mb8/group/model-boxcar/smoothed-4mm'
% each second-level SPM.mat can correspond to one contrast only!
% so now you need to run this separately for each contrast
matlabbatch{1}.spm.stats.factorial_design.dir = {out_path};
% now that these are contrast niftis they don't have a time dimension
% so each image only gets mentioned once, with ',1' appended. just do that in here
% hell yeah append() is vectorized for cell array input
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = append(paths_con, ',1');
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
% per above, each group contrast contains only one contrast already signed
% so we only need dummy contrasts for the positive and negative (flipped) versions of said
consess = cell(1, 2);  
consess{1}.tcon.name = contrast_name;
consess{1}.tcon.weights = 1;
consess{1}.tcon.sessrep = 'none';
consess{2}.tcon.name = append(contrast_name, '_inverse');
consess{2}.tcon.weights = -1;
consess{2}.tcon.sessrep = 'none';
% as usual, now append into the matlabbatch object
matlabbatch{3}.spm.stats.con.consess = consess;
matlabbatch{3}.spm.stats.con.delete = 1;

%% execute order 66
spm_jobman('run',matlabbatch);
