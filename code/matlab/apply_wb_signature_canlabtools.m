%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% VARIABLES THAT MUST BE DEFINED BEFORE THE SCRIPT IS SOURCED
% THIS SCRIPT IS DESIGNED TO RUN ACROSS A SINGLE TASK
% aggregating over subject-level 1st-level SPM beta or contrast statmaps
% each of these should come in from targets calling script

% paths_nifti: derived from target. 
% cell array of paths to the same beta/contrast nifti but different subjects
% out_path: for the table of cosine similarities with each signature (col) by subject (row). one char array for one output file.

%% GET SUBJECT NUMBERS FROM PATHS_NIFTI
% so it will be robust to potential differences in fold structure/order
subj_nums = zeros(length(paths_nifti), 1);
for i=1:length(paths_nifti)
    path_parts = strsplit(paths_nifti{i}, "/");
    % THIS IS HARD CODED TO WORK ON SPM FIRST LEVEL BETAS SAVED IN OUR SPECIFIC FOLDER STRUCTURE
    subj_parts = strsplit(path_parts{10}, "-");
    subj_num = str2double(subj_parts{2});
    subj_nums(i) = subj_num;
end

%% LOAD NIFTIS AS CANLABTOOLS FMRI_DATA

contrasts = fmri_data(paths_nifti);

%% LOAD CANLABTOOLS WHOLE BRAIN SIGNATURE SET
% this is hard-coded under the hood to pull from wherever Neuroimaging_Pattern_Masks lives
wb_signature = load_image_set('multiaversive');

%% APPLY WHOLE BRAIN SIGNATURE TO NIFTIS IN QUESTION
% It's so insane that canlabtools puts this similarity calculation thing into apply_mask
% But... whatever, man.
% this returns a matrix of subjects (niftis) on the rows x signature file on the col
signature_similarity = apply_mask(contrasts, wb_signature, 'pattern_expression', 'cosine_similarity');
signature_similarity = array2table(signature_similarity, VariableNames=cellstr(wb_signature.image_names));
signature_similarity.subj_num = subj_nums;

%% WRITE OUT
writetable(signature_similarity, out_path)
