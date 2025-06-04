classify_controlled <- function (cons_level1.5, outcome_var) {
  out <- cons_level1.5 %>% 
    rename(outcome = {{outcome_var}}) %>% 
    # x_prefix = "voxel" works for both bold and encoding (assuming encoding is always by space)
    fit_model_xval(x_prefix = "X",
                   y_prefix = "outcome",
                   parsnip_model = set_engine(parsnip::pls(mode = "regression", num_comp = 5), engine = "plsr", method = "simpls"),
                   additional_recipe_steps = \(x) step_dummy(x, all_of("outcome"), role = "outcome", one_hot = TRUE, skip = TRUE),
                   rm_x = TRUE,
                   return_longer = FALSE) %>% 
    select(preds) %>% 
    mutate(preds = map(preds, \(x) tidy_plsda_preds(x),
                       .progress = "cleaning up class predictions"),
           subj_num = map_dbl(preds, \(x) unique(x$subj_num)),
           preds = map(preds, \(x) select(x, -subj_num)))
}
