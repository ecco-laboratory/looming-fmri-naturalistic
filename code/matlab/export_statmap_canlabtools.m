%% library thangs
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore/CanlabCore/'))
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks/'))
addpath('/home/data/eccolab/Code/GitHub/spm12/')

%% variables that must be specified before the script is called
% region: a cell of 1 or more ROI names in NeuroimagingPatternMasks syntax. will all be treated as a union ROI if multiple are specified.
% zs: column vector for the stat values to be patched into the nifti
% must have come from canlabtools cause it's going back into canlabtools
% ALREADY THRESHOLDED IN THE TARGET SCRIPT. ALPHA IS THERE
% out_path: the filename for the _single_ statmap nifti output

%% load in and mask sample fmri data (to get the volInfo)
% set constant sample nifti using one subject's first-level beta image
% should be fine because the non-dynamic scan parameters of controlled and naturalistic are the same
% the logic below only works (?) when we patch the data into an image that was originally 3D, 
% so we can't do this with the first vol of a timeseries nifti
sample_data = fmri_data('/home/data/eccolab/SPLaT_fMRI/ignore/models/task-naturalistic/acq-mb8/sub-0001/model-boxcar/smoothed-4mm/beta_0001.nii');
% apparently need to patch the 0-value voxels so they don't get extra-masked out by apply_mask
% fine bc we just need the mask for shape and size anyway. the actual values will come from data
sample_data.dat(sample_data.dat == 0) = -99;
%% construct the statmap by ROI and combine into single object

% statistic_image is an object constructor from canlabtools
statmap = statistic_image;

% if ROI specified, pre-mask data so that the voxel dim lines up with the incoming data. else leave it as whole-brain
if exist('region', 'var') == 1
    % this script currently only writes out one ROI at a time to nifti bc idc
    masked_data = apply_mask(sample_data, select_this_atlas_subset(region));
else
    masked_data = sample_data;
end

% I guess, just to be safe, construct these the way Phil had them
% where the data info comes from an actual masked nifti of fmri data
statmap.volInfo = masked_data.volInfo;
statmap.removed_voxels = masked_data.removed_voxels;
fprintf('number of voxels remaining after mask: %d\n', sum(~masked_data.removed_voxels))
fprintf('number of voxels expected from input data: %d\n', length(zs))
% REMEMBER!! NEEDS TO BE COLUMN VECTOR
statmap.dat = zs;

%% prepare statistic_image object for writing out to nifti

statmap.fullpath = out_path;
disp('preparing to write statmap nifti')
write(statmap, 'overwrite');
