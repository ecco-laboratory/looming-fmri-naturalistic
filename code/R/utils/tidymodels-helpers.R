# targets-safe tidymodels wrapper functions ----
# to fit cross-validated models small enough not to be run in matlab

# assumes that all predictor cols and all outcome cols share a respective prefix. easy enough
make_recipe <- function (in_data, x_prefix, y_prefix) {
  in_recipe <- recipe(in_data) %>% 
    update_role(starts_with(x_prefix), new_role = "predictor") %>% 
    update_role(starts_with(y_prefix), new_role = "outcome") %>% 
    update_role(subj_num, new_role = "ID")
  
  return (in_recipe)
}

make_workflow <- function (in_data, x_prefix, y_prefix, parsnip_model) {
  in_recipe <- make_recipe(in_data, x_prefix, y_prefix)
  
  out_workflow <- workflow() %>% 
    add_model(parsnip_model) %>% 
    add_recipe(in_recipe)
  
  return (out_workflow)
}

# wrapper around predict() that does some additional post-processing of the predictions
get_preds_one_fold <- function (model, split, x_prefix, y_prefix, rm_x) {
  
  pred_only <- model %>% 
    predict(new_data = testing(split))
  
  if (ncol(pred_only) == 1) {
    # so that it matches the same naming pattern of the multi-outcome preds
    new_col_name <- paste0(".pred_", y_prefix)
    pred_only %<>%
      rename(!!new_col_name := ".pred")
  } 
  
  preds <- testing(split)
  
  if (rm_x) {
    preds %<>%
      dplyr::select(-starts_with(x_prefix))
  }
  preds %<>% 
    rename_with(\(x) paste0(".obs_", x), starts_with(y_prefix)) %>% 
    bind_cols(pred_only) %>% 
    pivot_longer(cols = c(starts_with(".obs"), 
                          starts_with(".pred")), 
                 names_to = c(".value", NA, y_prefix), 
                 names_sep = "_")
  
  return (preds)
}

fit_pls_xval <- function (in_data, 
                      x_prefix, 
                      y_prefix = "rating",
                      num_comp = 1,
                      rm_x = FALSE) {
  
  parsnip_model <- parsnip::pls(mode = "regression", num_comp = num_comp) %>% 
    set_engine(engine = "plsr",
               # same algorithm as matlab plsregress
               method = "simpls")
  
  this_workflow <- make_workflow(in_data, x_prefix = x_prefix, y_prefix = y_prefix, parsnip_model = parsnip_model)
  
  out_xval <- in_data %>% 
    group_vfold_cv(group = subj_num) %>% 
    mutate(model = purrr::map(splits, \(x) fit(this_workflow, data = training(x)), .progress = "fitting models"),
           preds = purrr::map2(model, splits, \(x, y) get_preds_one_fold(x, y, x_prefix = x_prefix, y_prefix = y_prefix, rm_x = rm_x), .progress = "getting predictions"))
  
  return (out_xval)
}