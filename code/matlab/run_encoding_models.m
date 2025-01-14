%% load libraries
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));

%% define relevant variables
fmri_dir = '/archival/projects/SPLaT/data/fmri/nifti';
derivatives_dir = fullfile(fmri_dir, 'derivatives/fmriprep-23.1.4');
activations_dir = '/home/data/eccolab/SPLaT_fMRI/ignore/outputs';
subjects = readtable(fullfile(fmri_dir, 'participants.tsv'), FileType="delimitedtext");
subjects = subjects(strcmp(subjects.group, 'use'), :);
task_defaults = readstruct('/home/data/eccolab/SPLaT_fMRI/task_defaults.json');
% just for the naturalistic task rn
task_defaults = task_defaults(2);

models ={'alexnet', 'flynet', 'flyalexnet'};
regions = {{'Bstem_SC'},{'Amy'}};

%% call the internal script if you must. RIGHT NOW THIS DOESN'T RUN!
if false
    trs_to_use = (1:task_defaults.n_trs_kept) + floorDiv(task_defaults.disdaq_duration, task_defaults.tr_duration);
    
    paths_nifti = cell(height(subjects), 1);
    paths_confounds = cell(height(subjects), 1);
    paths_activations = repmat({paths_nifti}, length(models), 1);
    
    % construct the inputs expected by the script bc not calling from targets
    % here, construct the stuff that's consistent by subject/encoding model
    for s=1:height(subjects)
        for r=1:task_defaults.n_runs
            
            paths_nifti{s}{r} = fullfile(derivatives_dir, subjects.participant_id{s}, 'func', sprintf('smoothed_4mm_%s_task-naturalistic_run-%02d_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii', subjects.participant_id{s}, r));
            paths_confounds{s}{r}= fullfile(derivatives_dir, subjects.participant_id{s}, 'func', sprintf('%s_task-naturalistic_run-%02d_desc-confounds_trimmed.txt', subjects.participant_id{s}, r));
            for m=1:length(models)   
                paths_activations{m}{s}{r}=fullfile(activations_dir, sprintf('task.naturalistic_%s_run.%02d_acts-%s.csv', strrep(subjects.participant_id{s}, '-', '.'), r, models{m}));  
            end
        end
    end
    
    tr_duration = task_defaults.tr_duration;
    
    fit_pls_canlabtools;
end

%% read in pred obs corrs that have been previously saved out
pred_obs_corrs = cell(length(models), 1);

for i=1:length(models)
    pred_obs_corrs{i} = cell(length(regions), 1);
    for j=1:length(regions)
        switch regions{j}{1}
            case 'Bstem_SC'
            region_fname = 'sc';
            case 'Amy'
            region_fname = 'amyg';
        end
        pred_obs_corrs{i}{j} = readmatrix(fullfile(activations_dir, sprintf('naturalistic_perf.%s_%s.csv', models{i}, region_fname)));
    end
end

%% load fMRI template and masks for visualization

template_path = fullfile(derivatives_dir, 'sub-0001/func/smoothed_4mm_sub-0001_task-naturalistic_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii');
template = fmri_data(template_path);

masks = cell(length(regions), 1);
for i=1:length(regions)
    masks{i} = select_this_atlas_subset(regions(i));
end

%% figures
for r=1:length(regions)
    
    for m=1:length(models)
        plot_data = apply_mask(template, masks{r});
        plot_data.dat = pred_obs_corrs{m}{r}';
        plot_data.removed_images = 0;
        tmap = ttest(plot_data, .05, 'FDR');
        orthviews(tmap)
    end
    
    
    plot_data.dat = pred_obs_corrs{2}{r}' - pred_obs_corrs{1}{r}';
    plot_data.removed_images = 0;
    tmap = ttest(plot_data,.05,'FDR');
    orthviews(tmap)
    
end

% 2025-01-14 MT: Please be warned. Nothing below this point works because it now depends on stuff that is not written out by the wrapper script!
% If you really want the whole brain connectivity it must be re-loaded
%% plot betas as a function of distance from a physical plane in the brain?
for m=1:length(models)
    for r=1:length(regions)
        plot_data = apply_mask(template, masks{r});
        plot_data.dat = pred_obs_corrs{m}{r}';
        plot_data.removed_images = 0;
        tmap = ttest(plot_data, .05, 'FDR');
        
        % set to false because neither yhat and all subjs whole-brain BOLD are getting written out of fit_pls_canlabtools atm. deal with it
        if false
            for s=1:n_subjs
                r_map(m,:,s) = corr(mean(squeeze(yhat_bymodel(m,subj_indices==s,:)),2),bold_allsubjs(subj_indices==s,:));
            end
        end
        
        % this section appears to be for calculating distance from an xy plane
        xyz = plot_data.volInfo.xyzlist*template.volInfo.mat(1:3,1:3);
        xyz(plot_data.removed_voxels,:)=[];
        for i=1:3
            xyz(:,i)=xyz(:,i) + plot_data.volInfo.mat(i,4);
        end
        
        xyz(:,1) =abs(xyz(:,1));
        
        distance_from_boundary=zeros(1,height(xyz));
        for i=1:height(xyz)
            distance_from_boundary(i) = pdist([[3 -30];xyz(i,1:2)]);
            
        end
        
        subplot(2,1,m);hold all
        title(models{m})
        for s=1:size(plot_data.dat,2)
            b(s,m,:) = glmfit(zscore(distance_from_boundary),template.dat(:,s));
            scatter(distance_from_boundary,plot_data.dat(:,s),'.')
        end
        
        xlabel 'Distance from [x = ±3, y = -30] (mm)'
        ylabel 'Beta coefficient'
    end
end

%% connectivity with average predicted response
figure;
template = fmri_data(template_path);
template.dat = squeeze(r_map(1,:,:));
tmap_conn = ttest(template,.01,'FDR');
montage(tmap_conn)
%% connectivity with average predicted response
figure;
template = fmri_data(template_path);
template.dat = squeeze(r_map(2,:,:));
tmap_conn = ttest(template,.01,'FDR');
montage(tmap_conn)
%% connectivity with average predicted response
figure;
template = fmri_data(template_path);
template.dat = squeeze(r_map(2,:,:)+r_map(1,:,:))/2;
tmap_conn = ttest(template,.01,'FDR');
montage(tmap_conn)
%% difference in connectivity; looming vs alexnet
figure;
template = fmri_data(template_path);
template.dat = squeeze(r_map(2,:,:)-r_map(1,:,:));
tmap_conn = ttest(template,.01,'FDR');
montage(tmap_conn)
