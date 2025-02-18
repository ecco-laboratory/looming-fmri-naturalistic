%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% variables that must be specified before the script is called
% model_path: full path to the level 1 SPM.mat (thus the output will be for 1 subj)
% out_folder: the folder stem for each output file containing data from one roi

%% extract beta condition names from the SPM.mat
% remember the SPM object is always saved under the varname `SPM`
load(model_path);
[model_folder, ~, ~] = fileparts(model_path);
beta_names = SPM.xX.name;
% the cell wrap syntax is necessary bc Vbeta is a struct array with length == n_betas
% these come out with only the file name, not the path
beta_files = {SPM.Vbeta.fname};
beta_files = fullfile(model_folder, beta_files);

% TODO: is it possible to get the nifti directly out of the SPM.mat?? check the docs
% this is stupid but all the actual condition names have more characters
beta_files = beta_files(strlength(beta_names) > 14);
beta_names = beta_names(strlength(beta_names) > 14);
% then strip 'Sn(N)' prefixes from the condition names (otherwise counterbalancing means condition names won't match)
beta_names = extractAfter(beta_names, 6);
beta_names = erase(beta_names, '*bf(1)');
% then reorder the files into condition alphabetical order instead of scan presentation order
% so they will be consistent across Ss
[beta_names, sort_order] = sort(beta_names);
beta_files = beta_files(sort_order);

%% actually do the processing
% this is made to be called from targets, so it just naked does the stuff
% as opposed to defining a matlab function to be called from within matlab
% select_atlas_subset expects a cell array
% can have more than one ROI and it will be the union of all the ROIs

% for a whole-brain-esque analysis, ROIs are now hard-set in this script!
% atm, go through all glasser parcels plus amygdala and superior colliculus
% TODO: annoyingly, the canlab2018 atlas doesn't breakdown amygdala by hemisphere. 
% only by nuclei but each nucleus ROI has both hemispheres?
% newer canlab atlases do break it up by hemisphere but no longer by nucleus (prob fine)
% consider asking about cloning a new copy of Neuroimaging_Pattern_Masks
atlas = load_atlas('canlab2018');
% bear in mind that this separates each of the Ctx parcels by hemisphere
rois = atlas.labels(startsWith(atlas.labels, 'Ctx'));
% but these are bilateral (SC has the option to split by hemi but by putting only the prefix canlabtools will group)
rois = [rois {'Amygdala'} {'Bstem_SC'}];

% fmri_data can be called on a cell array of 3D niftis
data=fmri_data(beta_files);

for i=1:length(rois)
    disp(['Current ROI: ' rois{i}])
    % critically, when it applies a mask to a stack of niftis
    % it will keep any voxel that is within the mask for ANY image
    % so we don't have to handle when some subjects have zero voxels in the mask
    masked_dat = apply_mask(data, select_this_atlas_subset(rois(i), atlas));
    % transpose to get file on the row and voxel on the column
    % apply the beta/condition names as row names
    % this format is better for PLS bc the voxels need to be 'variables'
    % however, when these are fed into R to generate RDMs, R cor() needs condition in the columns
    % if you want to use cor() to generate RDMs of the betas, 
    % you will need to transpose to get voxels in columns
    betas = array2table(masked_dat.dat', RowNames=beta_names);
    % the parcel name will be stored in the file name
    writetable(betas, fullfile(out_folder, sprintf('betas_%s.csv', rois{i})), WriteRowNames=true);
end
