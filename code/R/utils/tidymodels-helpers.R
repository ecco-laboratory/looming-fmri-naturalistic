# targets-safe tidymodels wrapper functions ----
# to fit cross-validated models small enough not to be run in matlab

# assumes that all predictor cols and all outcome cols share a respective prefix. easy enough
make_recipe <- function (in_data, x_prefix, y_prefix, additional_steps = NULL) {
  in_recipe <- recipe(in_data) %>% 
    update_role(starts_with(x_prefix), new_role = "predictor") %>% 
    update_role(starts_with(y_prefix), new_role = "outcome") %>% 
    update_role(subj_num, new_role = "ID") # %>% 
    # step_normalize(all_outcomes(), skip = TRUE)
  
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
                            num_comp = 1,
                            rm_x = FALSE) {
  
  parsnip_model <- parsnip::pls(mode = "regression", num_comp = num_comp) %>% 
    set_engine(engine = "plsr",
               # same algorithm as matlab plsregress
               method = "simpls")
  
  this_workflow <- make_workflow(in_data, x_prefix = x_prefix, y_prefix = y_prefix, parsnip_model = parsnip_model) %>% 
    fit(data = in_data)
  
  out_data <- get_preds_one_fold(this_workflow, in_data, x_prefix = x_prefix, y_prefix = y_prefix, rm_x = rm_x)
  
  return (out_data)
}

fit_model_2level <- function (in_data,
                              nest_vars,
                              x_prefix, 
                              y_var,
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
