%this script is an example of multivariate pathway identification in Monica's data from SPLaT;
for s = 1:15 %for each of 15 subjects
%% load an estimated model from SPM
load(['/home/data/eccolab/SPLaT_fMRI/ignore/models/task-controlled/acq-mb8/sub-' sprintf('%04d',s) '/model-boxcar/smoothed-4mm/SPM.mat']);

%% find 4D niis from fmriprep
files = dir(['/home/data/eccolab/SPLaT_fMRI/ignore/data/fmri/derivatives/fmriprep-23.1.4/sub-' sprintf('%04d',s) '/func/smoothed_4mm_sub-' sprintf('%04d',s) '_task-controlled_run-*_space-MNI152NLin2009cAsym*']);

%% loop over runs and concatenate into one large fmri data object
for f =1:length(files) 
    dat = fmri_data([files(f).folder filesep files(f).name]); 
    if f==1; DAT = dat; DAT.dat = DAT.dat(:,17:end); else; DAT.dat = [DAT.dat dat.dat(:,17:end)];
    end
end

%% get masks for IT, LIP, and SC
atl = load_atlas('canlab2018');
sc = select_atlas_subset(atl,'Bstem_SC');
glasser = load_atlas('Glasser');
LIP = select_atlas_subset(glasser,{'Ctx_LIP'});
IT = select_atlas_subset(glasser,{'Ctx_TE' , 'Ctx_VV'});

%% estimate the pathway models
[~,inds]=max(SPM.xX.X(:,end-4:end),[],2);
stats = model_brain_pathway(DAT, IT, LIP, sc, sc, 'Indices', inds);

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

latent_corrs(s,:)=mean(stats.latent_correlations);
end