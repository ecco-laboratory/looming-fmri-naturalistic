%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% variables that must be specified before the script is called
% rois: a cell array of 1 or more ROI names in NeuroimagingPatternMasks syntax
% can be a nested cell array to group multiple atlas ROI names into the same value in the output atlas
% out_path: the filename for the _single_ atlas-volume file

%% apply ze mask

atlas = load_atlas('canlab2018');

% this is made to be called from targets, so it just naked does the stuff
% as opposed to defining a matlab function to be called from within matlab
% select_atlas_subset expects a cell array
% can have more than one ROI and it will be the union of all the ROIs
% with multiple ROIs, it uses atlas-style encoding where each voxel within an ROI
% gets labeled with a positive integer
% fence post case to start the first ROI in the sequence
mask_all_rois = select_this_atlas_subset(rois{1}, atlas);
mask_all_rois.dat(mask_all_rois.dat ~= 0) = 1;
for i=2:length(rois)
    this_mask = select_this_atlas_subset(rois{i}, atlas);
    % re-label each ROI that is made up of a group of atlas ROIs 
    % to have the same number label
    this_mask.dat(this_mask.dat ~= 0) = i;
    % recursively add this ROI into the mega mask
    mask_all_rois = image_math(mask_all_rois, this_mask, 'add');
end

%% write ze nifti
% the built in canlabtools write() method writes to the fullpath field
mask_all_rois.fullpath = out_path;
write(mask_all_rois, 'overwrite');
