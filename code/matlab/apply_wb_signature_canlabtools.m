%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS A SINGLE TASK
% aggregating over subject-level 1st-level SPM beta or contrast statmaps
% each of these should come in from targets calling script

% paths_nifti: derived from target.
% length 1 cell array of single path to across-subjects statmap nifti, OR cell array of paths to the same beta/contrast nifti but different subjects.
% the calling function sets this as a cell array no matter the number of file paths.
% path_fmri_data: path to CSV input file with subjects (or whatever) on the rows and fmri_data whole brain voxel on the columns
% must have come from canlabtools cause it's going back into canlabtools
% DO NOT SET paths_nifti AND path_fmri_data AT THE SAME TIME. USE ONLY ONE!!! the calling script at least will error out if you set both
% image_set_name: char keyword for a set of signatures known to canlabtools load_image_set(). 
% the calling script currently defaults this to 'multiaversive' for Ceko 2022
% out_path: for the table of cosine similarities with each signature (col) by subject (row). one char array for one output file.

%% LOAD NIFTIS AS CANLABTOOLS FMRI_DATA, OR POP TABULAR DATA INTO FMRI_DATA OBJECT

if exist('paths_nifti', 'var') == 1
    data = fmri_data(paths_nifti);
elseif exist('path_fmri_data', 'var') == 1
    data = fmri_data('/home/data/eccolab/SPLaT_fMRI/ignore/models/task-naturalistic/acq-mb8/sub-0001/model-boxcar/smoothed-4mm/beta_0001.nii');
    wb_data = readmatrix(path_fmri_data);
    data.dat = wb_data';
end

%% LOAD CANLABTOOLS WHOLE BRAIN SIGNATURE SET
% this is hard-coded under the hood to pull from wherever Neuroimaging_Pattern_Masks lives
wb_signature = load_image_set(image_set_name);

%% APPLY WHOLE BRAIN SIGNATURE TO NIFTIS IN QUESTION
% this returns a matrix of subjects (niftis) on the rows x signature file on the col
signature_similarity = apply_mask(data, wb_signature, 'pattern_expression', 'ignore_missing', 'cosine_similarity');
signature_similarity = array2table(signature_similarity, VariableNames=cellstr(wb_signature.image_names));

%% IF VALID, GET SUBJECT NUMBERS FROM PATHS_NIFTI AND APPEND TO SIMILARITY VALUES

if exist('paths_nifti', 'var') == 1 & length(paths_nifti) > 1
    % so it will be robust to potential differences in fold structure/order
    subj_nums = zeros(length(paths_nifti), 1);
    for i=1:length(paths_nifti)
        path_parts = strsplit(paths_nifti{i}, "/");
        % THIS IS HARD CODED TO WORK ON SPM FIRST LEVEL BETAS SAVED IN OUR SPECIFIC FOLDER STRUCTURE
        subj_parts = strsplit(path_parts{10}, "-");
        subj_num = str2double(subj_parts{2});
        subj_nums(i) = subj_num;
    end
    
    signature_similarity.subj_num = subj_nums;
end

%% WRITE OUT
writetable(signature_similarity, out_path)
