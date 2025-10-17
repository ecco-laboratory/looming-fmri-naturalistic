function contrasts = calc_contrasts_level1_controlled(condition_names)
    contrasts = cell(1, 24);
    contrasts{1}.tcon.name = 'attend_animal';                                                                                                            
    contrasts{1}.tcon.weights = contains(condition_names, "animal") - contains(condition_names, "hemifield");
    contrasts{2}.tcon.name = 'dog';
    contrasts{2}.tcon.weights = contains(condition_names, "dog") - (contains(condition_names, "frog") + contains(condition_names, "spider"))/2;
    contrasts{3}.tcon.name = 'frog';
    contrasts{3}.tcon.weights = contains(condition_names, "frog") - (contains(condition_names, "dog") + contains(condition_names, "spider"))/2;
    contrasts{4}.tcon.name = 'spider';
    contrasts{4}.tcon.weights = contains(condition_names, "spider") - (contains(condition_names, "dog") + contains(condition_names, "frog"))/2;
    contrasts{5}.tcon.name = 'above';
    contrasts{5}.tcon.weights = contains(condition_names, "above") - contains(condition_names, "below");
    contrasts{6}.tcon.name = 'looming';
    contrasts{6}.tcon.weights = contains(condition_names, "looming") - contains(condition_names, "receding");
    contrasts{7}.tcon.name = 'looming_baseline';
    % the others get converted to double by being mathed
    % if no math, must manually convert
    contrasts{7}.tcon.weights = double(contains(condition_names, "looming"));
    % this block of 6 contrasts is to get level 1 maps averaged across above/below, for cross-analysis with the naturalistic task
    contrasts{8}.tcon.name = 'looming_dog_baseline';
    contrasts{8}.tcon.weights = double(contains(condition_names, "looming") & contains(condition_names, "dog"));
    contrasts{9}.tcon.name = 'looming_frog_baseline';
    contrasts{9}.tcon.weights = double(contains(condition_names, "looming") & contains(condition_names, "frog"));
    contrasts{10}.tcon.name = 'looming_spider_baseline';
    contrasts{10}.tcon.weights = double(contains(condition_names, "looming") & contains(condition_names, "spider"));
    contrasts{11}.tcon.name = 'receding_dog_baseline';
    contrasts{11}.tcon.weights = double(contains(condition_names, "receding") & contains(condition_names, "dog"));
    contrasts{12}.tcon.name = 'receding_frog_baseline';
    contrasts{12}.tcon.weights = double(contains(condition_names, "receding") & contains(condition_names, "frog"));
    contrasts{13}.tcon.name = 'receding_spider_baseline';
    contrasts{13}.tcon.weights = double(contains(condition_names, "receding") & contains(condition_names, "spider"));
    contrasts{14}.tcon.name = 'stimuli';
    contrasts{14}.tcon.weights = contains(condition_names, "looming") + contains(condition_names, "receding");
    contrasts{15}.tcon.name = 'ratings';
    contrasts{15}.tcon.weights = double(contains(condition_names, "ratings"));
    % Teehee: This uses the SPM session means to estimate voxelwise mean signal a la marsbar
    % to use for subsequent percent signal change calculations
    contrasts{16}.tcon.name = 'meansignal';
    contrasts{16}.tcon.weights = double(contains(condition_names, "constant"));
    
    %% baseline for all of the attend to object category trials
    contrasts{17}.tcon.name = 'attendObject_baseline';
    contrasts{17}.tcon.weights = double(contains(condition_names, "animal"));
    %% baseline for all of the attend to spatial location trials
    contrasts{18}.tcon.name = 'attendLocation_baseline';
    contrasts{18}.tcon.weights = double(contains(condition_names, "location"));

    %% 
    contrasts{19}.tcon.name = 'attendObject_dog_baseline';
    contrasts{19}.tcon.weights = double(contains(condition_names, "animal") & contains(condition_names, "dog"));
    contrasts{20}.tcon.name = 'attendObject_frog_baseline';
    contrasts{20}.tcon.weights = double(contains(condition_names, "animal") & contains(condition_names, "frog"));
    contrasts{21}.tcon.name = 'attendObject_spider_baseline';
    contrasts{21}.tcon.weights = double(contains(condition_names, "animal") & contains(condition_names, "spider"));

    %%
    contrasts{22}.tcon.name = 'attendLocation_dog_baseline';
    contrasts{22}.tcon.weights = double(contains(condition_names, "location") & contains(condition_names, "dog"));
    contrasts{23}.tcon.name = 'attendLocation_frog_baseline';
    contrasts{23}.tcon.weights = double(contains(condition_names, "location") & contains(condition_names, "frog"));
    contrasts{24}.tcon.name = 'attendLocation_spider_baseline';
    contrasts{24}.tcon.weights = double(contains(condition_names, "location") & contains(condition_names, "spider"));

end
