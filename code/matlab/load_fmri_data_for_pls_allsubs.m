function [fmri_data_allsubjs, subj_indices, subj_nums] = load_fmri_data_for_pls_allsubs(paths_fmri)
    n_subjs = length(paths_fmri);

    % fencepost for first subject to allow preallocation for all subjects
    bold = load(paths_fmri{1});
    % assume this is going to be true for every subject. same timeseries duration
    bold_height = height(bold.DATA);
    
    fmri_data_allsubjs = zeros(n_subjs*bold_height, width(bold.DATA));
    subj_indices = zeros(n_subjs*bold_height, 1);
    this_subj_indices = 1:bold_height;
    fmri_data_allsubjs(this_subj_indices, :) = bold.DATA;
    
    % pull the proper subject number from paths_masked and use THAT instead of loose fold order for subj indices
    % so it will be robust to potential differences in fold structure
    subj_nums = [];
    [~, file_fmri, ~] = fileparts(paths_fmri{1});
    this_subj_num = str2double(extractBetween(convertCharsToStrings(file_fmri), 5, 8));
    subj_indices(this_subj_indices, 1) = this_subj_num;
    subj_nums = [subj_nums, this_subj_num];
    
    % loop over the rest of the subjects
    for i=2:n_subjs
        clear bold
        bold = load(paths_fmri{i});
        % whatever. I doubt failure
        this_subj_indices = (1:bold_height) + (bold_height * (i-1));
        fmri_data_allsubjs(this_subj_indices, :) = bold.DATA;
        [~, file_fmri, ~] = fileparts(paths_fmri{i});
        this_subj_num = str2double(extractBetween(convertCharsToStrings(file_fmri), 5, 8));
        subj_indices(this_subj_indices, 1) = this_subj_num;
        subj_nums = [subj_nums, this_subj_num];
    end
end
