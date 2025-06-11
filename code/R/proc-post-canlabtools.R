# various targets-safe helper functions for extra modeling stuff on canlabtools outputs ----

## voxelwise statmap data operations ----

# set all sub-threshold voxel values to 0
threshold_tvals_pre_statmap <- function (tvals,                                        
                                         threshold_t = NULL,
                                         threshold_p = NULL,
                                         p_adjust_method = "BH") {
  # thresholds are absolute values
  # if both are provided, use threshold_p
  if (!is.null(threshold_p)) {
    pvals <- tvals %>% 
      # get two-tailed p-values
      map_dbl(\(x) map_pval_from_tval(x))
    
    if (!is.na(p_adjust_method)) pvals <- p.adjust(pvals, method = p_adjust_method, n = sum(!is.nan(pvals)))
    tvals[pvals > threshold_p] <- 0
    
  } else if (!is.null(threshold_t)) {
    tvals[abs(tvals) < threshold_t] <- 0
  }
  
  return (tvals)
}

# this function assumes voxels comes in on the COLUMN
# and held-out subject is on the ROW
# fun_compare must take 2 vector arguments, compare them vectorized/element-wise, and return 1 vector
# fun_compare is ignored if in_data_2 is not specified (and there is thus nothing to compare)
summarize_tvals_pre_statmap <- function (in_data_1, in_data_2 = NULL, fun_compare = magrittr::subtract) {
  
  in_data_1 <- as.list(in_data_1)
  
  if (!is.null(in_data_2)) {
    in_data_2 <- as.list(in_data_2)
    test_data <- map2(in_data_1, in_data_2, \(x, y) fun_compare(x, y),
                      .progress = "calculating comparison map over voxels")
  } else {
    test_data <- in_data_1
  }
  
  test_data %>% 
    # use list and naked map so we get the .progress bar. whole-brain can take a while bc so many voxels
    # get the cross-validated t-value for each voxel as the independent mean/SE over HELD-OUT SUBJECTS
    # leaving na.rm = FALSE for all of them so it only returns a connectivity t-val if all subjects had values
    map_dbl(\(x) mean(x) / (sd(x)/sqrt(length(x))), 
            .progress = "calculating voxel t-vals over xval folds") %>% 
    # because matlab doesn't have an NA type
    replace_na(NaN) %>% 
    # round so that sending these to matlab through the command line won't be too long
    round(digits = 3)
}

## modeling other stuff as a function of main fMRI encoding model preds ----

### testing naturalistic encoding model-predicted pattern on controlled task betas ----

calc_controlled_pattern_discrim <- function (patterns_controlled,
                                             path_activations_controlled,
                                             metadata_controlled,
                                             path_betas_encoding,
                                             betas_discrim_encoding) {
  activations <- read_activations_controlled(path_activations_controlled, metadata_controlled)
  
  betas_encoding <- read_betas_encoding(path_betas_encoding)
  
  patterns <- patterns_controlled %>% 
    calc_naturalistic_encoding_controlled_pred(activations, betas_encoding) %>% 
    left_join(betas_discrim_encoding %>% 
                rename(betas_discrim = coefs),
              by = "subj_num") %>% 
    mutate(discrim_pred = map2(pattern_pred, betas_discrim, \(x, y) {
      out <- c(x %*% y)
      names(out) <- colnames(y)
      return(out[c("outcome_dog", "outcome_frog", "outcome_spider")])
      }),
      animal_type_pred = map_chr(discrim_pred, \(x) str_split_i(names(which.max(x)), "_", 2)),
      across(starts_with("animal_type"), \(x) factor(x, levels = c("dog", "frog", "spider")))) %>% 
    select(subj_num, animal_type, direction, discrim_pred, animal_type_pred) %>% 
    unnest_wider(discrim_pred)
  
  return (patterns)
}

calc_controlled_perf <- function (pattern_expression_controlled, cor_method = "pearson") {
  pattern_expression_controlled %>% 
    select(-activations, -betas) %>% 
    unchop(cols = c(pattern, pattern_pred)) %>% 
    group_by(direction, subj_num) %>% 
    summarize(correlation = cor(pattern, pattern_pred, method = cor_method), .groups = "drop") %>% 
    pivot_wider(names_from = direction, values_from = correlation) %>% 
    mutate(diff = looming - receding) %>% 
    summarize(across(c(looming, receding, diff),
                     .fns = summary_stats_default),
              cohens.d = cohens.d(looming, receding))
}

calc_controlled_pattern_expression <- function (patterns_controlled,
                                                path_activations_controlled,
                                                metadata_controlled,
                                                path_betas_encoding) {
  
  activations <- read_activations_controlled(path_activations_controlled, metadata_controlled)
  
  betas_encoding <- read_betas_encoding(path_betas_encoding)
  
  patterns <- patterns_controlled %>% 
    calc_naturalistic_encoding_controlled_pred(activations, betas_encoding) %>% 
    mutate(pattern_expression = map2_dbl(pattern_pred, pattern, \(x, y) x %*% y),
           pred_obs_correlation = map2_dbl(pattern_pred, pattern, \(x, y) cor(x, y, method = "pearson")),
           roi_avg_pred = map_dbl(pattern_pred, mean))
  
  return (patterns)
}

calc_naturalistic_encoding_controlled_pred <- function (patterns_controlled,
                                                        activations_controlled,
                                                        betas_encoding) {
  out <- patterns_controlled %>% 
    # Matlab's PLS algorithm mean-centers X and Y column-wise before fitting
    # so do the same here to keep these within the expected scale
    # the scale is important because the pattern expression statistic is a product, and so sign-sensitive
    mutate(across(starts_with("X"), \(x) c(scale(x, scale = FALSE)))) %>% 
    nest(.by = c(subj_num, animal_type, direction), .key = "pattern") %>% 
    left_join(activations_controlled, by = c("animal_type", "direction")) %>% 
    left_join(betas_encoding, by = "subj_num") %>% 
    mutate(pattern = map(pattern, unlist),
           # colMeans on the product will mean over video frame in the rows
           pattern_pred = map2(activations, betas, \(x, y) colMeans(x %*% y)))
  
  return (out)
}

read_activations_controlled <- function (path_activations_controlled,
                                         metadata_controlled) {
  out <- read_csv(path_activations_controlled) %>% 
    right_join(metadata_controlled %>% 
                 select(filename, animal_type, direction),
               by = c("video" = "filename")) %>% 
    # average over videos within animal x direction condition
    # bc the controlled beta patterns are only at that resolution
    group_by(animal_type, direction, frame) %>% 
    summarize(across(where(is.numeric), mean), .groups = "drop") %>% 
    # from here on out you better not do anything to reorder the rows
    select(-frame) %>% 
    nest(.by = c(animal_type, direction), .key = "activations") %>% 
    mutate(activations = map(activations, as.matrix))
  
  return (out)
}

read_betas_encoding <- function (path_betas_encoding) {
  out <- read_csv(path_betas_encoding, col_names = FALSE) %>% 
    rename(subj_num = X1) %>% 
    nest(.by = subj_num, .key = "betas") %>% 
    # drop the intercept term
    mutate(betas = map(betas, \(x) slice_tail(x, n=nrow(x)-1)),
           betas = map(betas, as.matrix))
}

### classifying human-labeled stimulus category by encoding model-predicted pattern ----

compare_betas_encoding_category_selfreport <- function (betas_loom,
                                                        betas_object,
                                                        betas_selfreport,
                                                        n_bootstraps) {
  # get all sets of betas as tibbles with outcome on the column and voxel on the row
  betas_loom %<>%
    as_tibble() %>% 
    # in 2 class, they are direct opposite so only need one
    select(outcome_loom)
  
  betas_object %<>%
    as_tibble()
  
  betas_selfreport %<>% 
    unnest_longer(coefs) %>% 
    pivot_wider(names_from = rating_type, values_from = coefs) %>% 
    select(-coefs_id) %>% 
    rename_with(\(x) paste0("rating_", x), everything()) %>% 
    # so that the self report betas will all correlate positively with each other
    mutate(rating_unpleasantness = -rating_pleasantness) %>% 
    select(-rating_pleasantness)

  boots <- bind_cols(betas_loom, betas_object, betas_selfreport) %>% 
    # not ideal but bootstrap at the order of voxels because bootstrapping earlier 
    # at the level of re-fitting the second-order models is not really tractable
    bootstraps(times = n_bootstraps, apparent = TRUE) %>% 
    mutate(correlations = map(splits,
                              \(x) analysis(x) %>% 
                                cor() %>% 
                                as_tibble(rownames = "row") %>% 
                                pivot_longer(cols = -row, names_to = "col", values_to = "correlation"),
                              .progress = "bootstrapping cross-outcome beta cormats")) %>% 
    select(-splits)
  
  apparent <- boots %>% 
    filter(id == "Apparent") %>% 
    unnest(correlations) %>% 
    select(-id) %>% 
    rename(cor_true = correlation)
  
  out <- boots %>% 
    filter(id != "Apparent") %>% 
    unnest(correlations) %>% 
    group_by(row, col) %>% 
    summarize(across(correlation,
                     .fns = summary_stats_default,
                     .names = "cor_{.fn}"), 
              .groups = "drop") %>% 
    left_join(apparent, by = c("row", "col"))

  return (out)
}

calc_perm_pval_object_by_pattern <- function (preds, perms, grouping_cols = NULL) {
  
  grouping_cols <- enquo(grouping_cols)
  acc_true <- preds %>% 
    select(preds) %>% 
    unnest(preds) %>% 
    calc_metrics_nested_classprobs(grouping_cols = !!grouping_cols)
  
  # don't need this until a few lines down but do it before I quote the argument
  join_by_cols <- c(".metric", ".estimator")
  if (!quo_is_null(grouping_cols)) join_by_cols <- c(join_by_cols, quo_name(grouping_cols))
  
  out <- perms %>% 
    unnest(acc_perm) %>% 
    left_join(acc_true, by = join_by_cols, suffix = c("_perm", "_real")) %>% 
    group_by(.metric, pick(!!grouping_cols)) %>% 
    summarize(estimate_real = unique(.estimate_real),
              # need this to report permuted chance levels
              estimate_perm = median(.estimate_perm),
              pval = (sum(.estimate_perm > .estimate_real)+1)/(n()+1))
  
  return (out)
}

permute_acc_object_by_pattern <- function (preds, n_perms, outcome_categories = c("object", "looming", "object_looming"), acc_grouping_cols = NULL) {
  stopifnot(length(outcome_categories) == 1)
  acc_grouping_cols <- enquo(acc_grouping_cols)
  
  out <- preds %>% 
    # keep only cols required for computing accuracy
    select(subj_num, run_num, !!acc_grouping_cols, .obs, .preds, .pred) %>% 
    # nest blocks that will be kept together in permuting
    nest(obs = .obs, pred = c(!!acc_grouping_cols, .preds, .pred), .by = c(subj_num, run_num)) %>% 
    # the n_trs stuff is to shorten all "runs" to the shortest one so that the permuted TR vectors will be the same length
    # because omitting the fixation TRs causes the runs no longer to have the same amount of data within or between subject
    mutate(n_trs = map_int(obs, nrow)) %>% 
    group_by(subj_num) %>% 
    mutate(min_n_trs = min(n_trs)) %>% 
    ungroup() %>% 
    # shorten by dropping later TRs
    mutate(obs = map2(obs, min_n_trs, \(x, y) slice_head(x, n = y)),
           pred = map2(pred, min_n_trs, \(x, y) slice_head(x, n = y))) %>% 
    # next, by subject so that blocks will only be permuted within subject
    nest(.by = subj_num) %>%
    # IF THIS IS CALLED WITHIN A TAR_REP, MULTIPLY N_PERMS BY THE NUMBER OF BATCHES FOR THE TOTAL OVERALL NUMBER OF ITERATIONS
    mutate(data = map(data, \(x) permutations(x, obs, times = n_perms), 
                      .progress = "permuting")) %>% 
    # now get splits out of the nest and into a column
    unnest(data) %>% 
    # reconstruct across-subjects permuted datasets. requires pulling full data out of the splits objects
    mutate(perm_data = map(splits, \(x) analysis(x) %>% unnest(cols = c(obs, pred)), 
                           .progress = "unnesting within subject")) %>% 
    select(-splits) %>% 
    unnest(perm_data) %>% 
    nest(.by = id)
  
  out %<>% 
    # now compute permuted accuracy & AUROC metrics
    mutate(acc_perm = map(data, \(x) calc_metrics_nested_classprobs(x, grouping_cols = !!acc_grouping_cols), 
                          .progress = "calculating perm-wise acc")) %>% 
    # drop the actual permuted data for space saving
    select(-data)
  
  return (out)
}

fit_object_by_pattern <- function (path_pattern_allsubs, 
                                   events_allsubs, 
                                   n_trs_kept_per_run, 
                                   path_pattern_allsubs_2 = NULL,
                                   pattern_type = c("bold", "encoding"), 
                                   outcome_categories = c("obj", "loom", "obj.loom"),
                                   n_pls_comp = 5,
                                   xval = TRUE) {
  # just to save memory since these aren't used for this analysis
  events_allsubs %<>%
    select(-starts_with("rating"))
  
  stopifnot(length(pattern_type) == 1)
  
  if (pattern_type == "bold") {
    
    timecourses <- path_pattern_allsubs %>% 
      load_bold_mat_allsubs() %>% 
      select(-fold_num)
  
  } else if (pattern_type == "encoding") {
    
    timecourses <- load_encoding_pred_allsubs(path_pattern_allsubs)
    
    if (!is.null(path_pattern_allsubs_2)) {
      timecourses %<>%
        left_join(load_encoding_pred_allsubs(path_pattern_allsubs_2),
                  by = c("subj_num", "tr_num"),
                  suffix = c(".1", ".2"))
    }
  }
  
  out <- timecourses %>% 
    left_join(events_allsubs, by = c("subj_num", "tr_num")) %>% 
    # reconstruct run number here for later block permutation by run
    mutate(run_num = (tr_num-1) %/% n_trs_kept_per_run) %>% 
    # PHIL HAD SET UP HIS VERSION OF THE ANALYSIS TO EXCLUDE FIXATION AND FOOD TIMEPOINTS
    # since food never looms and is very different from the animals
    filter(!is.na(animal_type), !(animal_type %in% c("fixation", "food")))
  
  stopifnot(length(outcome_categories) == 1)
  
  if (outcome_categories == "obj") {
    out %<>%
      rename(outcome = animal_type)
  } else if (outcome_categories == "loom") {
    out %<>%
      mutate(outcome = if_else(has_loom == 1, "loom", "no.loom"))
  } else if (outcome_categories == "obj.loom") {
    out %<>%
      unite("outcome", animal_type, has_loom, sep = ".")
  }
  
  if (!xval) {
    out %<>% 
      # fit_pls_single returns a list object
      fit_pls_single(x_prefix = "voxel",
                     y_prefix = "outcome",
                     num_comp = n_pls_comp,
                     additional_recipe_steps = \(x) step_dummy(x, all_of("outcome"), role = "outcome", one_hot = TRUE, skip = TRUE),
                     rm_x = TRUE,
                     return_longer = FALSE) %>% 
      pluck("model") %>% 
      # now JUST the coefs
      extract_pls_workflow_coefs()
    
    # if no xval, we aren't looking at any held-out predictions. just the coefs
    return (out)
  }
  
  out %<>%
    # x_prefix = "voxel" works for both bold and encoding (assuming encoding is always by space)
    fit_model_xval(x_prefix = "voxel",
                   y_prefix = "outcome",
                   parsnip_model = set_engine(parsnip::pls(mode = "regression", num_comp = n_pls_comp), engine = "plsr", method = "simpls"),
                   additional_recipe_steps = \(x) step_dummy(x, all_of("outcome"), role = "outcome", one_hot = TRUE, skip = TRUE),
                   rm_x = TRUE,
                   return_longer = FALSE) %>% 
    # 2025-05-05: We want to save the betas out now too for subsequent second-order correlation analysis with the encoding.selfreport models
    # the object is too large to save out if we keep all the model fits
    mutate(coefs = map(model, \(x) extract_pls_workflow_coefs(x),
                       .progress = "extracting coefs")) %>% 
    select(coefs, preds) %>% 
    mutate(preds = map(preds, \(x) tidy_plsda_preds(x),
                       .progress = "cleaning up class predictions"),
           subj_num = map_dbl(preds, \(x) unique(x$subj_num)),
           preds = map(preds, \(x) select(x, -subj_num)))
  
  return (out)
}

get_ratings_by_encoding_time <- function (path_pred_allsubs, events_allsubs) {
  
  pred_encoding <- load_encoding_pred_allsubs(path_pred_allsubs)
  
  out <- pred_encoding %>% 
    # ??? SHOULD WE INDEPENDENTLY Z-SCORE EACH OF THE PRED-BOLD COLUMNS BEFORE AVERAGING TOGETHER?
    mutate(across(starts_with("voxel"), \(x) c(scale(x)))) %>% 
    pivot_longer(cols = starts_with("voxel"),
                 names_to = "voxel_num",
                 values_to = "bold_pred") %>% 
    group_by(subj_num, tr_num) %>% 
    summarize(bold_pred_mean = mean(bold_pred), 
              .groups = "drop") %>% 
    left_join(events_allsubs, by = c("subj_num", "tr_num")) %>% 
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

get_ratings_by_encoding_space <- function (path_pred_allsubs, events_allsubs, path_pred_allsubs_2 = NULL) {
  
  pred_encoding <- load_encoding_pred_allsubs(path_pred_allsubs)
  
  if (!is.null(path_pred_allsubs_2)) {
    pred_encoding %<>%
      left_join(load_encoding_pred_allsubs(path_pred_allsubs_2),
                by = c("subj_num", "tr_num"),
                suffix = c(".1", ".2"))
  }
  
  out <- pred_encoding %>% 
    left_join(events_allsubs, by = c("subj_num", "tr_num")) %>% 
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
  
  encoding_preds_1 <- load_encoding_pred_allsubs(path_encoding_preds_1)
  
  encoding_preds_2 <- load_encoding_pred_allsubs(path_encoding_preds_2)
  
  coefs_xval <- encoding_preds_1 %>% 
    full_join(encoding_preds_2, by = c("subj_num", "tr_num"), suffix = paste0("_", encoding_names_str)) %>% 
    full_join(bold, by = c("subj_num", "tr_num")) %>% 
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
    rename(subj_num = X1) %>% 
    rename_with(\(x) x %>% 
                  str_sub(start = 2L) %>% 
                  as.integer() %>% 
                  subtract(1L) %>% 
                  paste0("voxel.", .), 
                starts_with("X")) %>% 
    group_by(subj_num) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup()
  
  return (out)
}
