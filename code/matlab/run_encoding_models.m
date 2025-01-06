ss = [1:10 12:15 17 19 ];
confound_files = dir('/archival/projects/SPLaT/data/fmri/nifti/derivatives/fmriprep-23.1.4/**/func/sub-0*naturalistic*_trimmed.txt');
files = dir('/archival/projects/SPLaT/data/fmri/nifti/derivatives/fmriprep-23.1.4/**/func/smoothed_4mm_sub-0*naturalistic*.nii');
models ={'alexnet','flynet'};
regions = {'Bstem_SC','Amy'};
trs_to_use = 17:1005;

    for rr=1:length(regions);
for m=1:length(models);

        cc=1;

        for s=1:16;
            for r=1:3;
                paths_nifti{s}{r}=[files(cc).folder filesep files(cc).name];
                paths_confounds{s}{r}= [confound_files(cc).folder filesep confound_files(cc).name];
                paths_activations{s}{r}=sprintf(['/home/data/eccolab/SPLaT_fMRI/ignore/outputs/task.naturalistic_sub.%04d_run.%02d_acts-' models{m} '.csv'],ss(s),r);

                cc=cc+1;
            end;
        end
        tr_duration = .492;
        region = regions{rr};    paths_activations(5)=[];
        paths_confounds(5)=[];
        paths_nifti(5)=[];

        fit_pls_canlabtools;

        activations_bymodel{rr,m} = activations;
        yhat_bymodel{rr,m} =yhat;
        pred_obs_corr_bymodel{rr,m} =pred_obs_corr;
    end
end

%% combine alexnet and flynet models


for k=1:n_subjs
    train_idx = subj_indices~=k;
    test_idx = ~train_idx;
    [~,~,~,~,beta_cv] = plsregress(zscore([squeeze(activations_bymodel{1}(train_idx,:)) squeeze(activations_bymodel{2}(train_idx,:)) ]), bold_masked_allsubjs(train_idx,:), 100);

    % ones() prepends an intercept column to the selected chunk of conv_features to be used for predicting against the betas
    yhat_combined(test_idx,:)=[ones(sum(test_idx),1) zscore([squeeze(activations_bymodel{1}(test_idx,:)) squeeze(activations_bymodel{2}(test_idx,:)) ])]*beta_cv;

    % has to be like this bc corr(X, Y) correlates each pair of columns (here, voxels)
    % but we only care about correlating each voxel's real data to its own predicted data
    % so diag() keeps only each voxel to itself
    pred_obs_corr_combined(k,:)=diag(corr(yhat_combined(test_idx,:), bold_masked_allsubjs(test_idx,:)));

end


%%

template_path = '/home/data/eccolab/SPLaT_fMRI/ignore/data/fmri/derivatives/fmriprep-23.1.4/sub-0001/func/smoothed_4mm_sub-0001_task-naturalistic_run-01_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii';

template = fmri_data(template_path);
template=apply_mask(template,mask);

for m=1:2
    template.dat = squeeze(pred_obs_corr_bymodel(m,:,:))';
    template.removed_images = 0;
    tmap = ttest(template,.05,'FDR');

    for s=1:n_subjs
        r_map(m,:,s) = corr(mean(squeeze(yhat_bymodel(m,subj_indices==s,:)),2),bold_allsubjs(subj_indices==s,:));
    end


    xyz = template.volInfo.xyzlist*template.volInfo.mat(1:3,1:3);
    xyz(template.removed_voxels,:)=[];
    for i=1:3

        xyz(:,i)=xyz(:,i)+template.volInfo.mat(i,4);

    end

    xyz(:,1) =abs(xyz(:,1));

    distance_from_boundary=zeros(1,height(xyz));
    for i=1:height(xyz)
        distance_from_boundary(i) =  pdist([[3 -30];xyz(i,1:2)]);

    end

    subplot(2,1,m);hold all
    title(models{m})
    for s=1:size(template.dat,2)
        b(s,m,:) = glmfit(zscore(distance_from_boundary),template.dat(:,s));
        scatter(distance_from_boundary,template.dat(:,s),'.')
    end

    xlabel 'Distance from [x = ±3, y = -30] (mm)'
    ylabel 'Beta coefficient'

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


%% figures

template = fmri_data(template_path);
template=apply_mask(template,mask);

template.dat = squeeze(pred_obs_corr_bymodel(1,:,:))';
template.removed_images = 0;
tmap_alexnet = ttest(template,.05,'FDR');
orthviews(tmap_alexnet)

template.dat = squeeze(pred_obs_corr_bymodel(2,:,:))';
template.removed_images = 0;
tmap_flynet = ttest(template,.05,'FDR');
orthviews(tmap_flynet)


template.dat = squeeze(pred_obs_corr_bymodel(2,:,:))' -squeeze(pred_obs_corr_bymodel(1,:,:))';
template.removed_images = 0;
tmap_flynet = ttest(template,.05,'FDR');
orthviews(tmap_flynet)
