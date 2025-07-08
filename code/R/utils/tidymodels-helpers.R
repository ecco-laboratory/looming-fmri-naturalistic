# targets-safe tidymodels wrapper functions ----
# to fit cross-validated models small enough not to be run in matlab

# assumes that all predictor cols and all outcome cols share a respective prefix. easy enough
make_recipe <- function (in_data, x_prefix, y_prefix, additional_steps = NULL) {
  in_recipe <- recipe(in_data) %>% 
    update_role(starts_with(x_prefix), new_role = "predictor") %>% 
    update_role(starts_with(y_prefix), new_role = "outcome") %>% 
    update_role(subj_num, new_role = "ID")
  
  if (!is.null(additional_steps)) {
    in_recipe %<>%
      additional_steps()
  }
  
  return (in_recipe)
}

make_workflow <- function (in_data, x_prefix, y_prefix, parsnip_model, additional_recipe_steps = NULL) {
  in_recipe <- make_recipe(in_data, x_prefix, y_prefix, additional_steps = additional_recipe_steps)
  
  out_workflow <- workflow() %>% 
    add_model(parsnip_model) %>% 
    add_recipe(in_recipe)
  
  return (out_workflow)
}

# always extracts from the solution using the maximal number of components
extract_pls_workflow_coefs <- function (workflow_fit) {
  model_obj <- extract_fit_engine(workflow_fit)
  ncomp <- pluck(model_obj, "ncomp")
  coefficients <- pluck(model_obj, "coefficients")
  
  return(drop(coefficients[ , , ncomp]))
}

# wrapper around predict() that does some additional post-processing of the predictions
get_preds_one_fold <- function (model, data_or_split, x_prefix, y_prefix, rm_x, return_longer = TRUE) {
  if ("rsplit" %in% class(data_or_split)) {
    pred_only <- model %>% 
      predict(new_data = testing(data_or_split))
    
    preds <- testing(data_or_split)
  } else if ("data.frame" %in% class(data_or_split)) {
    pred_only <- model %>% 
      predict(new_data = data_or_split)
    
    preds <- data_or_split
  }
  
  if (ncol(pred_only) == 1) {
    # so that it matches the same naming pattern of the multi-outcome preds
    new_col_name <- paste0(".pred_", y_prefix)
    current_col_name <- names(pred_only)
    pred_only %<>%
      rename(!!new_col_name := !!current_col_name)
  } 
  
  if (rm_x) {
    preds %<>%
      dplyr::select(-starts_with(x_prefix))
  }
  preds %<>% 
    rename_with(\(x) paste0(".obs_", x), starts_with(y_prefix)) %>% 
    bind_cols(pred_only)
  
  if (return_longer) {
    preds %<>% 
      pivot_longer(cols = c(starts_with(".obs"), 
                            starts_with(".pred")), 
                   names_to = c(".value", NA, y_prefix), 
                   names_sep = "_")
  }
  
  return (preds)
}

tidy_plsda_preds <- function (preds, apply_softmax = FALSE) {
  out <- preds %>% 
    rename(.obs = .obs_outcome) %>% 
    # at this point, .obs is just the original categorical column but .pred has one columns for each one-hot outcome level
    # bundle .obs and .pred together to make it easier to extract the max for each one as the class label
    nest(.preds = starts_with(".pred")) %>% 
    mutate(.preds = map(.preds, \(y) unlist(y)),
           # keep both the numeric PLS-DA "class predictions" and the label for the max class prediction
           # the former for ROC curves and the latter for straight-up accuracy
           .pred = map_chr(.preds, \(y) names(which.max(y))))
  
  if (apply_softmax) {
    out %<>%
      mutate(.preds = map(.preds, \(y) softmax(y)))
  }
  
  out %<>%
    mutate(# once the max class prediction has been pulled,
      # keep only the first class prob value for 2-class scenarios
      # because tidymodels only takes one class prob for 2-class auc
      .preds = map(.preds, \(y) if (length(y) == 2) return (y[1]) else return (y)),
      # clean up class labels that were pulled from predictor col names
      .pred = str_split_i(.pred, "_", 3),
      # and finally push everything back to factor for tidymodels accuracy etc
      .obs = factor(.obs),
      # patch all categories back in in case one is never predicted
      .pred = factor(.pred, levels = levels(.obs)))
}

fit_model_xval <- function (in_data, 
                            x_prefix, 
                            y_prefix = "rating",
                            additional_recipe_steps = NULL,
                            parsnip_model = set_engine(parsnip::pls(mode = "regression", num_comp = 1, engine = "plsr"), method = "simpls"),
                            rm_x = FALSE,
                            return_longer = TRUE) {
  
  this_workflow <- make_workflow(in_data, x_prefix = x_prefix, y_prefix = y_prefix, parsnip_model = parsnip_model, additional_recipe_steps = additional_recipe_steps)
  
  out_xval <- in_data %>% 
    group_vfold_cv(group = subj_num) %>% 
    mutate(model = purrr::map(splits, \(x) fit(this_workflow, data = training(x)), 
                              .progress = "fitting models"),
           preds = purrr::map2(model, splits, \(x, y) get_preds_one_fold(x, y, x_prefix = x_prefix, y_prefix = y_prefix, rm_x = rm_x, return_longer = return_longer), 
                               .progress = "getting predictions"))
  
  return (out_xval)
}

fit_pls_single <- function (in_data, 
                            x_prefix, 
                            y_prefix = "rating",
                            additional_recipe_steps = NULL,
                            num_comp = 1,
                            rm_x = FALSE,
                            return_longer = TRUE) {
  
  parsnip_model <- parsnip::pls(mode = "regression", num_comp = num_comp) %>% 
    set_engine(engine = "plsr",
               # same algorithm as matlab plsregress
               method = "simpls")
  
  this_workflow <- make_workflow(in_data, x_prefix = x_prefix, y_prefix = y_prefix, parsnip_model = parsnip_model, additional_recipe_steps = additional_recipe_steps) %>% 
    fit(data = in_data)
  
  out_data <- get_preds_one_fold(this_workflow, in_data, x_prefix = x_prefix, y_prefix = y_prefix, rm_x = rm_x, return_longer = return_longer)
  
  return (list(model = this_workflow, preds = out_data))
}

fit_model_2level <- function (in_data,
                              nest_vars,
                              x_prefix, 
                              y_var,
                              additional_recipe_steps = NULL,
                              parsnip_model,
                              n_folds_within_nest = 2,
                              rm_x = FALSE) {
  
  this_recipe <- in_data %>% 
    head() %>% 
    # because they won't be present in the nested data
    select(-any_of(nest_vars)) %>% 
    recipe() %>% 
    update_role(starts_with(x_prefix), new_role = "predictor") %>% 
    update_role(all_of(y_var), new_role = "outcome")
  
  if (!is.null(additional_recipe_steps)) {
    this_recipe %<>%
      additional_recipe_steps()
  }
  
  this_workflow <- workflow() %>% 
    add_model(parsnip_model) %>% 
    add_recipe(this_recipe)
  
  out <- in_data %>% 
    nest(.by = all_of(nest_vars)) %>% 
    mutate(cv = purrr::map(data, \(x) vfold_cv(x, v = n_folds_within_nest))) %>% 
    unnest(cv) %>% 
    mutate(model = purrr::map(splits, \(x) fit(this_workflow, data = training(x)), .progress = "fitting models"),
           preds = purrr::map2(model, splits, \(x, y) get_preds_one_fold(x, y, x_prefix = x_prefix, y_prefix = y_var, rm_x = rm_x), .progress = "getting predictions"))
  
  return (out)
}

calc_xval_perf <- function (df_preds, grouping_cols = NULL, classprob_prefix = ".pred_outcome", df_perms = NULL) {
  grouping_cols <- enquo(grouping_cols)
  join_by_cols <- ".metric"
  if (!quo_is_null(grouping_cols)) join_by_cols <- c(join_by_cols, quo_name(grouping_cols))
  
  acc_true <- df_preds %>% 
    select(subj_num, preds) %>% 
    # calculate accuracy within held-out subject
    mutate(acc = map(preds, \(x) calc_metrics_nested_classprobs(x,
                                                                grouping_cols = !!grouping_cols,
                                                                classprob_prefix = classprob_prefix))) %>% 
    select(-preds) %>% 
    unnest(acc) %>% 
    group_by(.metric, pick(!!grouping_cols)) %>% 
    summarize(across(.estimate, list(mean = mean, se = \(x) sd(x)/sqrt(length(x)))), .groups = "drop")
  
  if (!is.null(df_perms)) {
    out <- df_perms %>% 
      unnest(acc_perm) %>% 
      rename(.estimate_perm = .estimate) %>% 
      left_join(acc_true, by = join_by_cols) %>% 
      group_by(.metric, pick(!!grouping_cols)) %>% 
      # pval needs to be calculated first before the columns get overwritten with their summary values
      summarize(pval = (sum(.estimate_perm > .estimate_mean)+1)/(n()+1),
                across(c(.estimate_mean, .estimate_se), unique),
                # need this to report permuted chance levels
                .estimate_perm = median(.estimate_perm)) %>% 
      select(!!grouping_cols, .estimate_mean, .estimate_se, .estimate_perm, pval)
      
    
    return (out)
  } else {
    return (acc_true)
  }
}

permute_xval_classification <- function (df_preds, n_perms, grouping_cols = NULL, classprob_prefix = ".pred_outcome") {
  grouping_cols <- enquo(grouping_cols)
  
  # expects by-subj nested xval output to permute within subject
  out <- df_preds %>% 
    # IF THIS IS CALLED WITHIN A TAR_REP, MULTIPLY N_PERMS BY THE NUMBER OF BATCHES FOR THE TOTAL OVERALL NUMBER OF ITERATIONS
    mutate(splits_nested = map(preds, \(x) permutations(x, .obs, times = n_perms), 
                               .progress = "permuting")) %>% 
    select(-preds) %>% 
    # now get splits out of the nest and into a column
    unnest(splits_nested) %>% 
    # reconstruct across-subjects permuted datasets. requires pulling full data out of the splits objects
    mutate(perm_preds = map(splits, \(x) analysis(x), 
                            .progress = "unnesting within subject")) %>% 
    select(-splits) %>% 
    unnest(perm_preds) %>% 
    # this puts them back into "datasets" by permutation iteration where each iteration has all subjects
    nest(.by = id)
  
  out %<>% 
    # now compute permuted accuracy & AUROC metrics
    mutate(acc_perm = map(data, \(x) calc_metrics_nested_classprobs(x, 
                                                                    grouping_cols = !!grouping_cols,
                                                                    classprob_prefix = classprob_prefix), 
                          .progress = "calculating perm-wise acc")) %>% 
    # drop the actual permuted data for space saving
    select(-data)
  
  return (out)
}

calc_metrics_nested_classprobs <- function (preds, grouping_cols = NULL, classprob_prefix = ".pred_outcome") {
  grouping_cols <- enquo(grouping_cols)
  
  preds %>% 
    unnest_wider(col = .preds) %>% 
    group_by(pick(!!grouping_cols)) %>% 
    metrics(truth = .obs, estimate = .pred, starts_with(classprob_prefix))
}
