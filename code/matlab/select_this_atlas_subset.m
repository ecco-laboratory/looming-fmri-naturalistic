function mask = select_this_atlas_subset(region, canlab_atlas)
    
    addpath('/home/data/eccolab/Code/GitHub/spm12'); % per spm docs, do not genpath it
    addpath(genpath('/home/data/eccolab/Code/GitHub/CanlabCore'));
    addpath(genpath('/home/data/eccolab/Code/GitHub/Neuroimaging_Pattern_Masks'));
    % this atlas is not implemented through CanlabCore/Neuroimaging_Pattern_Masks 
    % so we're just gonna have to literally read in some niftis from this folder
    BrainstemNavigator_path = '/home/data/shared/BrainstemNavigator/0.9/2a.BrainstemNucleiAtlas_MNI/labels_thresholded_binary_0.35';
    
    % default atlas if none specified
    switch nargin
        case 1
            canlab_atlas = load_atlas('canlab2018');
        case 2
            if ischar(canlab_atlas)
                canlab_atlas = load_atlas(canlab_atlas);
            end
    end
    
    % load the mask for this ROI
    % for superior colliculus, manually use BrainstemNavigator ROI instead
    if any(strcmp(region,'Bstem_SC'))
        mask = fmri_data(fullfile(BrainstemNavigator_path, 'SC_l.nii'));
        mask_r = fmri_data(fullfile(BrainstemNavigator_path, 'SC_r.nii'));
        mask.dat = mask.dat + mask_r.dat;
        % and then explicitly exclude PAG from that SC ROI
        pag = load_atlas('Kragel2019PAG');
        pag = resample_space(pag,mask);
        mask.dat(pag.dat>0) = 0;
        clear mask_r pag
    else
        mask = select_atlas_subset(canlab_atlas, region);
    end
end
