%this script is an example of multivariate pathway identification in Monica's data from SPLaT;
%% add toolboxes to path
addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));


%% get masks for IT, LIP, and SC
atl = load_atlas('canlab2018');
% MT TODO: Feed the ROIs in as arguments from targets
sc = select_atlas_subset(atl,'Bstem_SC');
glasser = load_atlas('Glasser');
LIP = select_atlas_subset(glasser,{'Ctx_LIP'});
IT = select_atlas_subset(glasser,{'Ctx_TE' , 'Ctx_VV'});
pag = load_atlas('Kragel2019PAG');
vta = select_atlas_subset(load_atlas('CIT168'),{'VTA'});
MD = select_atlas_subset(atl,{'MD'});
FEF = select_atlas_subset(atl,{'FEF'});

TR = .492; %specify tr for filtering
%% loop over subjects and do analysis
for s = [1:10 12:15  ] %for each of 15 subjects; 1:10
    %% load an estimated model from SPM
    load(['/home/data/eccolab/SPLaT_fMRI/ignore/models/task-controlled/acq-mb8/sub-' sprintf('%04d',s) '/model-boxcar/smoothed-4mm/SPM.mat']);

    %% find 4D niis from fmriprep
    files = dir(['/home/data/eccolab/SPLaT_fMRI/ignore/data/fmri/derivatives/fmriprep-23.1.4/sub-' sprintf('%04d',s) '/func/smoothed_4mm_sub-' sprintf('%04d',s) '_task-controlled_run-*_space-MNI152NLin2009cAsym*']);
    confound_files = dir(['/home/data/eccolab/SPLaT_fMRI/ignore/data/fmri/derivatives/fmriprep-23.1.4/sub-' sprintf('%04d',s) '/func/sub-' sprintf('%04d',s) '_task-controlled_run-*_trimmed.txt']);

    %% loop over runs and concatenate into one large fmri data object
    for f =1:length(files)
        dat = fmri_data([files(f).folder filesep files(f).name]);
        % MT TODO: Implement a clearer way of discarding the first n volumes from the 8 second wait period
        if f==1;
            DAT = dat; DAT.dat = DAT.dat(:,17:end);
        else
            DAT.dat = [DAT.dat dat.dat(:,17:end)];
        end
    end


   
    % [rmssd, rmssd_outlier_regressor_matrix] = rmssd_movie(DAT,'showmovie',false,'nodisplay');

    % then read in and row-run-stack fmriprep motion regressors
    confounds = []; session_means =[];
    for j=1:length(confound_files)
        confounds = [confounds; readmatrix([confound_files(j).folder filesep confound_files(j).name])];
        session_means = [session_means; j*ones(height(readmatrix([confound_files(j).folder filesep confound_files(j).name])),1)];
    end


    preprocessed_dat = canlab_connectivity_preproc(DAT,'hpf', 1/128,TR, 'no_plots');




    %% estimate the pathway models
    % This is currently written with within-subject leave-one-run-out cross-validation
    
    stats = model_brain_pathway(preprocessed_dat, IT, LIP, sc, sc, 'Indices', session_means);
    stats_mb = model_brain_pathway(preprocessed_dat, pag, vta, sc, sc, 'Indices', session_means);
    stats_fef = model_brain_pathway(preprocessed_dat, MD, FEF, sc, sc, 'Indices', session_means);

    % %% create an object for visualization and look at the patterns in a source region
    % temp = stats.source_two_obj;
    % temp.dat = stats.Z_pathway_four;
    % orthviews(temp);
    %
    % %% create an object for visualization and look at the patterns in a target region
    % temp = stats.target_two_obj;
    % temp.dat = stats.V_pathway_four;
    % orthviews(temp);

    Z_lip{s}=stats.Z_pathway_four;
    Z_IT{s}=stats.Z_pathway_one;

    V_lip{s}=stats.V_pathway_four;
    V_IT{s}=stats.V_pathway_one;

    Z_VTA{s}=stats_mb.Z_pathway_four;
    Z_PAG{s}=stats_mb.Z_pathway_one;

    V_VTA{s}=stats_mb.V_pathway_four;
    V_PAG{s}=stats_mb.V_pathway_one;

     Z_md{s}=stats_fef.Z_pathway_one;
    Z_fef{s}=stats_fef.Z_pathway_four;

    V_md{s}=stats_fef.V_pathway_one;
    V_fef{s}=stats_fef.V_pathway_four;

    latent_corrs(s,:)=mean(stats.latent_correlations);
    latent_corrs_mb(s,:)=mean(stats_mb.latent_correlations);
    latent_corrs_fef(s,:)=mean(stats_fef.latent_correlations);

    b_md_sc(s,:)=glmfit(SPM.xX.X,stats_fef.latent_timeseries_pathway1(:,1));
    b_fef_sc(s,:)=glmfit(SPM.xX.X,stats_fef.latent_timeseries_pathway2(:,2));

    b_it_sc(s,:)=glmfit(SPM.xX.X,stats.latent_timeseries_pathway1(:,1));
    b_lip_sc(s,:)=glmfit(SPM.xX.X,stats.latent_timeseries_pathway2(:,2));

    b_vta_sc(s,:) =glmfit(SPM.xX.X,stats_mb.latent_timeseries_pathway1(:,1));
    b_pag_sc(s,:) =glmfit(SPM.xX.X,stats_mb.latent_timeseries_pathway2(:,2));



    X = [stats_fef.latent_timeseries_pathway1(:,1) stats_fef.latent_timeseries_pathway2(:,2) stats.latent_timeseries_pathway1(:,1) stats.latent_timeseries_pathway2(:,2) stats_mb.latent_timeseries_pathway1(:,1) stats_mb.latent_timeseries_pathway2(:,2)];

    wh_regressors = [ 1 1 0 0 0 0 ;...
        0 0 1 1 0 0;...
        0 0 0 0 1 1; ...
        1 1 1 1 1 1];

    for m=1:height(wh_regressors)

        for ii = 1:5
            train = session_means~=ii;
            test = ~train;
            [~,~,~,~,b_pls_second] = plsregress(zscore(X(train,wh_regressors(m,:)==1)),stats.target_one_obj.dat(:,train)');
            yhat(m,test,:)=[ones(length(find(test)),1) zscore(X(test,wh_regressors(m,:)==1))]*b_pls_second;
        end

        xval_prediction_outcome_correlation(s,m,:)=diag(corr(squeeze(yhat(m,:,:)),stats.target_one_obj.dat(:,:)'));
        x=squeeze(yhat(m,:,:))-stats.target_one_obj.dat'; x= x(:); %all voxels, all timepoints

        loglik(s,m) = -sum(log(normpdf(x,0,std(x))));
        bic(s,m) = -2*loglik(s,m)+length(find(wh_regressors(m,:)==1))*log(height(X)*width(X));
    end
end

%% compare performance of different pathway models

X=(squeeze(mean(xval_prediction_outcome_correlation,3)));
X(all(X'==0),:)=[];
barplot_columns(X,'dolines')
set(gca,'XTickLabel',{'FEF, MD','IT, LIP','PAG, VTA','Combined'})
ylabel("Voxel-wise Performance")
xlabel 'Pathways'


%%  create an object for visualization and look at the patterns in a source region
    
    % combined vs individual
    temp = stats_mb.target_two_obj;%stats.source_two_obj;
    temp.dat = squeeze(xval_prediction_outcome_correlation(:,4,:))'-squeeze(xval_prediction_outcome_correlation(:,1,:))'/3 -squeeze(xval_prediction_outcome_correlation(:,2,:))'/3 -squeeze(xval_prediction_outcome_correlation(:,3,:))'/3;
    temp.dat(:,all(temp.dat==0))=[];
    orthviews(threshold(ttest(temp),.01,'FDR'));


    % FEF/MD 
        temp = stats_mb.target_two_obj;%stats.source_two_obj;
    temp.dat = squeeze(xval_prediction_outcome_correlation(:,1,:))';
    temp.dat(:,all(temp.dat==0))=[];
     orthviews(threshold(ttest(temp),.01,'FDR'));


%LIP/IT 
       temp = stats_mb.target_two_obj;%stats.source_two_obj;
    temp.dat = squeeze(xval_prediction_outcome_correlation(:,2,:))';
    temp.dat(:,all(temp.dat==0))=[];
     orthviews(threshold(ttest(temp),.01,'FDR'));


    %PAG/VTA 
       temp = stats_mb.target_two_obj;%stats.source_two_obj;
    temp.dat = squeeze(xval_prediction_outcome_correlation(:,3,:))';
    temp.dat(:,all(temp.dat==0))=[];
     orthviews(threshold(ttest(temp),.01,'FDR'));


