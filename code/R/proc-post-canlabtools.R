# various targets-safe helper functions for extra modeling stuff on canlabtools outputs ----

# this function assumes voxels comes in on the COLUMN
# and held-out subject is on the ROW
summarize_tvals_pre_statmap <- function (in_data) {
  in_data %>% 
    # get the cross-validated t-value for each voxel as the independent mean/SE over HELD-OUT SUBJECTS
    summarize(across(everything(), 
                     \(x) mean(x) / (sd(x)/sqrt(length(x))) )) %>% 
    # convert the 1-row dataframe to a vector. the voxel values should stay in the same order
    as.numeric() %>% 
    # round so that sending these to matlab through the command line won't be too long
    round(digits = 3)
}

get_ratings_by_encoding_time <- function (path_pred_allsubs, events_allsubs) {
  
  pred_encoding <- load_encoding_pred_allsubs(path_pred_allsubs) %>% 
    group_by(fold_num) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup()
  
  out <- pred_encoding %>% 
    # ??? SHOULD WE INDEPENDENTLY Z-SCORE EACH OF THE PRED-BOLD COLUMNS BEFORE AVERAGING TOGETHER?
    mutate(across(starts_with("voxel"), \(x) c(scale(x)))) %>% 
    pivot_longer(cols = starts_with("voxel"),
                 names_to = "voxel_num",
                 values_to = "bold_pred") %>% 
    group_by(fold_num, tr_num) %>% 
    summarize(bold_pred_mean = mean(bold_pred), 
              .groups = "drop") %>% 
    left_join(events_allsubs, by = c("fold_num", "tr_num")) %>% 
    filter(!is.na(video_id)) %>% 
    group_by(subj_num, video_id) %>% 
    mutate(across(c(starts_with("rating"), animal_type, has_loom), \(x) x[1]),
           # to count off the TRs within each video endspike
           tr_num_temp = 1:n()) %>% 
    ungroup() %>% 
    select(-tr_num) %>% 
    # pivot time onto the columns to get however many endspike columns of bold pred per stimulus
    pivot_wider(names_from = tr_num_temp,
                values_from = bold_pred_mean,
                names_prefix = "tr_")
  
  return (out)
}

get_ratings_by_encoding_space <- function (path_pred_allsubs, events_allsubs) {
  
  pred_encoding <- load_encoding_pred_allsubs(path_pred_allsubs) %>% 
    group_by(fold_num) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup()
  
  out <- pred_encoding %>% 
    left_join(events_allsubs, by = c("fold_num", "tr_num")) %>% 
    filter(!is.na(video_id)) %>% 
    group_by(subj_num, video_id) %>% 
    mutate(across(c(starts_with("rating"), animal_type, has_loom), \(x) x[1])) %>% 
    # many fewer lines of code to keep voxel and average over TR because voxel was already on the column
    group_by(subj_num, fold_num, video_id, pick(starts_with("rating")), animal_type, has_loom) %>% 
    summarize(across(starts_with("voxel"), mean), .groups = "drop")
  
  return (out)
}

# taking in all subjs bold, and all subjs encoding preds from 2 candidate models to be partial-correlated with the bold against each other 
calc_pcor_encoding_on_bold <- function (path_bold_mat_all_subs, 
                                        path_encoding_preds_1, 
                                        path_encoding_preds_2,
                                        encoding_name_1,
                                        encoding_name_2) {
  
  encoding_names_str <- c(encoding_name_1, encoding_name_2)
  encoding_name_1 <- ensym(encoding_name_1)
  encoding_name_2 <- ensym(encoding_name_2)
  
  bold <- load_bold_mat_allsubs(path_bold_mat_all_subs) %>% 
    # ahead of joining it onto the voxel preds for each encoding type
    rename_with(\(x) paste0(x, "_obs"), starts_with("voxel"))
  
  encoding_preds_1 <- load_encoding_pred_allsubs(path_encoding_preds_1) %>% 
    group_by(fold_num) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup()
  
  encoding_preds_2 <- load_encoding_pred_allsubs(path_encoding_preds_2) %>% 
    group_by(fold_num) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup()
  
  coefs_xval <- encoding_preds_1 %>% 
    full_join(encoding_preds_2, by = c("fold_num", "tr_num"), suffix = paste0("_", encoding_names_str)) %>% 
    full_join(bold, by = c("fold_num", "tr_num")) %>% 
    pivot_longer(cols = starts_with("voxel"),
                 names_to = c("voxel_num", ".value"),
                 names_sep = "_") %>% 
    mutate(voxel_num = as.integer(str_sub(voxel_num, start = 7L))) %>% 
    # nest by HELD-OUT SUBJECT. currently treats voxels x timepoints as exchangeable within each single-subject lm
    nest(preds = -subj_num) %>% 
    # this is the workhorse
    # interpret the betas as partial correlations because all the Xs and Ys are z-scored
    # by a quick test outside, it's not an exact match, but it's suuuuper close. on the order of .005 correlation units different.
    mutate(coefs = map(preds, \(x) lm(scale(obs) ~ scale(!!encoding_name_1) + scale(!!encoding_name_2), data = x) %>% 
                         broom::tidy())) %>% 
    select(-preds) %>% 
    unnest(coefs) %>% 
    filter(term != "(Intercept)") %>% 
    # strip out the "scale()" from the predictor names
    mutate(term = str_sub(term, start = 7L, end = -2L))
  
  return (coefs_xval)
}

## functions to read in and preprocess matlab-created file targets ----
# expects the bold.masked.ROI_all.subs object, which is a list of paths to the individual subject masked timecourse mat files
load_bold_mat_allsubs <- function (path_bold_mat_allsubs) {
  out <- path_bold_mat_allsubs %>% 
    map(\(x) x %>% 
          R.matlab::readMat() %>% 
          pluck("DATA") %>% 
          as_tibble(.name_repair = "unique")) %>% 
    bind_rows(.id = "target_name") %>% 
    mutate(subj_num = as.integer(str_sub(target_name, start = -4L))) %>% 
    select(-target_name) %>% 
    rename_with(\(x) str_replace(x, "...", "voxel."), starts_with("...")) %>% 
    nest(bold = -subj_num) %>% 
    mutate(fold_num = 1:n(),
           bold = map(bold, \(x) mutate(x, tr_num = 1:nrow(x)))) %>% 
    unnest(bold) %>% 
    select(subj_num, fold_num, tr_num, everything())
  
  return (out)
}

load_encoding_pred_allsubs <- function (path_pred_allsubs) {
  out <- path_pred_allsubs %>% 
    read_csv(col_names = FALSE) %>% 
    rename(fold_num = X1) %>% 
    rename_with(\(x) x %>% 
                  str_sub(start = 2L) %>% 
                  as.integer() %>% 
                  subtract(1L) %>% 
                  paste0("voxel.", .), 
                starts_with("X"))
  
  return (out)
}
