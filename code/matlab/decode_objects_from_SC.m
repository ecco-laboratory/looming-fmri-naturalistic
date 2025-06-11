
models = {'alexnet','flynet'};
do_bootstrap = false;
for m = 1:length(models)
    
    %% get some information about file structure and number of subjects
    t = readtable(['naturalistic_betas.' models{m} '.only_sc.csv']);
    
    yhat = table2array(readtable(['naturalistic_pred.' models{m} '.only_sc.csv'])); %alexnet
    subjects = unique(yhat(:,1));
    subjects = subjects(1:end-3); %exclude a few for now
    categories = {'Dog','Cat','Frog','Spider','Food'};
    %% read in and organize data
    crun = [];cX=[];cY=[];cS=[]; cY_loom=[];
    for s = 1:length(subjects)
        run=[];
        for r=1:3
            events = readtable(['task-naturalistic_sub-00' sprintf('%02d',subjects(s)) '_run-0' num2str(r) '_events.csv']);
                
                if r==1 
                    t = readtable(['task-naturalistic_sub-00' sprintf('%02d',subjects(s)) '_run-0' num2str(r) '_acts-' models{m} '.csv']);
                    if s==1
                        cevents = events;
                        revents = cevents;
                    else
                        revents = events;
                    end
                else
                    t = [t; readtable(['task-naturalistic_sub-00' sprintf('%02d',subjects(s)) '_run-0' num2str(r) '_acts-' models{m} '.csv'])];
                    
                    cevents = [cevents; events];
                    revents = [revents;events];
                end
                
            run = [run; r*ones(height(readtable(['task-naturalistic_sub-00' sprintf('%02d',subjects(s)) '_run-0' num2str(r) '_events.csv'])),1)];
                
        end
        
        X = yhat(yhat(:,1)==subjects(s),2:end);
        
        cond_vec = strncmp(revents.condition,'dog',3)+2*strncmp(revents.condition,'cat',3)+3*strncmp(revents.condition,'frog',3)+4*strncmp(revents.condition,'spider',3)+5*strncmp(revents.condition,'food',3);
        
        
        
        cond_vec_loom = contains(revents.condition,'loom0')+2*contains(revents.condition,'loom1');
        
        Y_loom = condf2indic(cond_vec_loom);
        Y_loom = Y_loom(Y_loom(:,1)==0,:);
        Y_loom = Y_loom(:,2:end);
        
        cY_loom=[cY_loom;Y_loom];
        
        
        Y = condf2indic(cond_vec);
        run = run(Y(:,1)==0);
        
        X = X(Y(:,1)==0,:);
        Y = Y(Y(:,1)==0,:);
        Y = Y(:,2:end);
        
        crun=[crun;run];
        cX=[cX;X];
        cY=[cY;Y];
        cS = [cS; s*ones(height(X),1)];
    end
    
    %% get rid of food stimuli
    cX = cX(cY(:,end)==0,:);
    cS = cS(cY(:,end)==0);
    cY = cY(cY(:,end)==0,:);
    cY = cY(:,1:4);
    %% cross validation across subjects
    %  kinds = floor((1:height(run))/140); %10 fold
    % kinds = crossvalind('k',height(X),5); %random 5-fold
    
    kinds = cS;
    clear Yh
    for k=1:max(kinds);
        train = kinds~=k;
        test = ~train;
        [~,~,~,~,b]=plsregress(cX(train,:),cY(train,:),40);
        Yh(test,:) = [ones(length(find(test)),1) cX(test,:)]*b;
        [~,gt{k}]=max(cY(test,:),[],2);
        gt{k}=categories(gt{k});
        yh{k}=Yh(test,:);
        
        
        
    end
    
    Y_pred{m}=Yh(:,1:4);
    
    if m==length(models)
        for k=1:max(kinds);
            train = kinds~=k;
            test = ~train;
            r_part(k,:,:,1)=partialcorr(Y_pred{2}(test,:),cY(test,:),Y_pred{1}(test,:));
            r_part(k,:,:,2)=partialcorr(Y_pred{1}(test,:),cY(test,:),Y_pred{2}(test,:));
        end
    end
    [~,pred_vals{m}]=max(Yh,[],2);
    [~,gt_inds{m}]=max(cY,[],2);
    % roc{m} = rocmetrics(categories(gt_inds{m}),Yh,categories);
    
    for i=1:length(yh); 
        yh{i}=yh{i}(:,1:4);
    end
    roc{m}=rocmetrics(gt,yh,categories(1:4));
    
    
    %% confusion matrix (row normalized)
    cm=confusionmat(gt_inds{m},pred_vals{m});
    
    cm = bsxfun(@rdivide,cm,sum(cm,2));
    figure;imagesc(cm)
    set(gca,'XTickLabel',{'Dog','Cat','Frog','Spider','Food'})
    set(gca,'YTickLabel',{'Dog','Cat','Frog','Spider','Food'})
    set(gca,'XTick',1:5,'YTick',1:5)
    
    
    %% bootstrap pls coefficients
    
    if do_bootstrap
        [~,~,~,~,b]=plsregress(cX,cY,40);
        bs_b = bootstrp(1000,@bootpls,cX,cY);
        bs_b = reshape(bs_b,1000,size(b,1),[]);
        bs_mean = squeeze(mean(bs_b));
        bs_se = squeeze(std(bs_b));
        bs_Z = bs_mean ./ bs_se;
        bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
        
        
        sample_data = fmri_data('/home/data/eccolab/SPLaT_fMRI/ignore/models/task-naturalistic/acq-mb8/sub-0001/model-boxcar/smoothed-4mm/beta_0001.nii');
        statmap = statistic_image;
        masked_data = apply_mask(sample_data, select_this_atlas_subset('Bstem_SC'));
        
        statmap.volInfo = masked_data.volInfo;
        statmap.removed_voxels = masked_data.removed_voxels;
        statmap.dat = bs_mean(2:end,:);
        statmap.p = bs_P(2:end,:);
        orthviews(threshold(statmap,.05,'FDR','k',3))
    end
    
end