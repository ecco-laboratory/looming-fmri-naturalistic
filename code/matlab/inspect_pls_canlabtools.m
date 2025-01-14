%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% read in previously estimated PLS performance
pred_obs_corr = readmatrix('/home/data/eccolab/SPLaT_fMRI/ignore/outputs/naturalistic_perf.flynet_sc.csv');
mask = select_atlas_subset(load_atlas('canlab2018'), {'Bstem_SC'});
template_path = '/home/data/eccolab/SPLaT_fMRI/ignore/data/fmri/derivatives/fmriprep-23.1.4/sub-0001/func/smoothed_4mm_sub-0001_task-naturalistic_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii';

%% figures

template = fmri_data(template_path);
template=apply_mask(template,mask);
template.dat = pred_obs_corr';
template.removed_images = 0;
tmap = ttest(template,.01,'FDR');

create_figure('cutaways'); axis off

surface_handles = [addbrain('brainstem') addbrain('sc') addbrain('pag') addbrain('rn') addbrain('VTA')];

tmap.surface('surface_handles', surface_handles, 'noverbose','clim',[5 12]);

drawnow, snapnow