## targets-optimized functions for generating plots from data targets ----

plot_selfreport_ratings <- function (beh_data) {
  n_subjs <- length(unique(beh_data$subj_num))
  
  preplot <- beh_data %>% 
    mutate(has_loom = if_else(has_loom == 1, "Looming", "No looming"),
           animal_type = fct_relevel(animal_type, "food", "dog", "cat")) %>% 
    pivot_longer(cols = starts_with("rating"), names_to = "construct", values_to = "rating", names_prefix = "rating_") %>% 
    mutate(construct = if_else(construct == "pleasantness", "valence (- to +)", construct),
           construct = fct_relevel(construct, "valence (- to +)"))
  
  preplot %>% 
    # so that stuff will be plotted by subject
    group_by(subj_num, animal_type, has_loom, construct) %>%
    summarize(rating = mean(rating), .groups = "drop") %>%
    ggplot(aes(x = animal_type, y = rating, color = has_loom)) +
    # individual subjects spaghetti
    geom_line(aes(group = interaction(subj_num, has_loom)), alpha = 0.1) +
    # summary lines
    geom_line(aes(group = has_loom), stat = "summary", fun = "mean") +
    geom_pointrange(stat = "summary", 
                    fun.data = "mean_se") +
    labs(x = "Object type",
         y = "Self-report rating",
         color = NULL,
         subtitle = glue::glue("N = {n_subjs} participants")) +
    facet_wrap(~ construct)
  
}

# takes a long df target of cross-validated encoding predictions, a la pred.encoding.xval_{ROI}
plot_sample_timecourse_encoding <- function (pred_encoding_xval, encoding_models = NULL, subj = 1) {
  # only 1 subject at a time!
  stopifnot(length(subj) == 1)
  
  preplot <- pred_encoding_xval %>% 
    select(encoding_type, subj_num, everything()) %>% 
    # just do it for subject 1 by default because this is for a schematic
    # but set it as an arg so you COULD set it later if you needed to
    filter(subj_num == subj)
  
  if (!is.null(encoding_models)) {
    encoding_names <- names(encoding_models)
    preplot %<>%
      filter(encoding_type %in% encoding_models) %>% 
      # to relabel and reorder them into the names and level order specified by the argument. for plot facet order
      mutate(encoding_type = fct_recode(encoding_type, !!!encoding_models),
             encoding_type = fct_relevel(encoding_type, !!encoding_names))
  }
  
  preplot %<>%
    group_by(encoding_type) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup() %>% 
    pivot_longer(cols = starts_with("voxel"), 
                 names_to = "voxel_num", 
                 values_to = "bold_pred", 
                 names_transform = list(voxel_num = \(x) as.integer(str_split_i(x, "\\.", 2L)) - 1L)) %>%
    # just to cut down on what's getting drawn on the plot
    filter(voxel_num %% 5 == 0) %>% 
    ggplot(aes(x = tr_num, y = bold_pred, color = voxel_num)) + 
    geom_line(aes(group = voxel_num), alpha = 0.5) + 
    scale_color_viridis_b() + 
    guides(color = "none") + 
    facet_grid(encoding_type ~ .,
               scales = "free_y") +
    labs(x = "TR (concatenated across runs)",
         y = "encoding-model-predicted BOLD")
}

# takes a named list of 2 across-run single-subject model activation targets
# bc the other schematics also facet by encoding model type
plot_sample_timecourse_activation <- function (targets_subj_activation) {
  encoding_levels <- names(targets_subj_activation)
  
  targets_subj_activation %>% 
    # need to read in the sets of activations separately and keep them listed for a little while bc they don't have the same numbers of units (aka columns)
    map(\(x) bind_rows(map(x, \(y) read_csv(y, col_names = FALSE)))) %>% 
    map(\(x) mutate(x, tr_num = 1:nrow(x)) %>% 
          # that's also why the pivot is inside this map
          pivot_longer(cols = -tr_num, 
                       names_to = "unit_num", 
                       values_to = "activation") %>% 
          nest(.by = unit_num) %>% 
          # this has the effect of hard-coding the same number of selected units across models. 
          # not necessary bc it's already pivoted but fine
          slice_sample(n = 50)) %>% 
    bind_rows(.id = "encoding_type") %>% 
    group_by(encoding_type) %>% 
    mutate(unit_num_arbitrary = 1:n()) %>% 
    ungroup() %>% 
    unnest(data) %>% 
    mutate(encoding_type = fct_relevel(encoding_type, !!encoding_levels)) %>% 
    ggplot(aes(x = tr_num, y = activation, color = unit_num_arbitrary)) + 
    geom_line(aes(group = unit_num), alpha = 0.5) + 
    scale_color_viridis_b() + 
    guides(color = "none") +
    facet_grid(encoding_type ~ .,
               scales = "free_y") +
    labs(x = "TR (concatenated across runs)",
         y = "neural network model activation")
}

# expects a SINGLE-SUBJECT bold.masked target. since it was going to get filtered for one subject anyway
plot_sample_timecourse_bold <- function (path_bold_masked_mat) {
  path_bold_masked_mat %>% 
    R.matlab::readMat() %>% 
    pluck("DATA") %>% 
    as_tibble(.name_repair = "unique") %>% 
    mutate(tr_num = 1:n()) %>% 
    pivot_longer(cols = starts_with("..."), 
                 names_to = "voxel_num", 
                 values_to = "bold", 
                 names_transform = list(voxel_num = \(x) as.integer(str_sub(x, start = 4L)))) %>% 
    # again to cut down the number of voxels getting plotted
    filter(voxel_num %% 5 == 0) %>% 
    ggplot(aes(x = tr_num, y = bold, color = voxel_num)) + 
    geom_line(aes(group = voxel_num), alpha = 0.5) + 
    scale_color_viridis_b() + 
    guides(color = "none") +
    labs(x = "TR (concatenated across runs)",
         y = "BOLD")
}

plot_auroc_8cat <- function (object_preds) {
  object_preds %>% 
    unnest_wider(col = .preds) %>% 
    roc_curve(truth = .obs, starts_with(".pred_outcome")) %>% 
    separate_wider_delim(cols = .level, delim = ".", 
                         names = c("animal_type", "has_loom"), 
                         cols_remove = FALSE) %>% 
    mutate(has_loom = if_else(has_loom == 1, "Looming", "No looming")) %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = has_loom)) + 
    geom_path() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") + 
    coord_fixed() + 
    facet_wrap(~fct_relevel(animal_type, "dog")) +
    labs(color = NULL)
}

plot_confusion_8cat <- function(object_preds, facet_var = NULL) {
  facet_var <- enquo(facet_var)

  plot_data <- object_preds %>% 
    count({{facet_var}}, .obs, .pred) %>% 
    complete({{facet_var}}, .obs, .pred, fill = list(n = 0)) %>%
    group_by(pick({{facet_var}}, .obs)) %>% 
    mutate(prob = n/sum(n)) %>% 
    ungroup() %>% 
    mutate(across(c(.obs, .pred), \(x) fct_recode(x,
                                                  "cat, no looming" = "cat.0",
                                                  "cat, looming" = "cat.1",
                                                  "dog, no looming" = "dog.0",
                                                  "dog, looming" = "dog.1",
                                                  "frog, no looming" = "frog.0",
                                                  "frog, looming" = "frog.1",
                                                  "spider, no looming" = "spider.0",
                                                  "spider, looming" = "spider.1")),
           across(c(.obs, .pred), \(x) fct_relevel(x,
                                                   "dog, no looming",
                                                   "dog, looming")))
  
  plot_out <- plot_data %>% 
    ggplot(aes(x = .pred, y = fct_rev(.obs), fill = prob)) + 
    geom_raster() + 
    geom_text(aes(label = round(prob, 2)), color = "white") + 
    scale_fill_viridis_c(limits = c(0, .8)) +
    guides(x = guide_axis(angle = 30)) +
    labs(x = "Predicted class", y = "True class", fill = "Classification\nprobability") +
    facet_wrap(vars(!!facet_var))
  
  return (plot_out)
}

plot_encoding_performance <- function (perf_combined, encoding_types = NULL, rois = c("sc", "amyg")) {
  n_subjs <- length(unique(perf_combined$fold_num))
  
  perf <- perf_combined %>% 
    filter(roi %in% rois)
  
  if (!is.null(encoding_types)) {
    perf %<>%
      filter(encoding_type %in% encoding_types) %>% 
      mutate(encoding_type = fct_relevel(encoding_type, !!encoding_types),
             encoding_type = fct_recode(encoding_type, !!!encoding_types),
             encoding_family = encoding_type)
  } else {
    perf %<>%
      mutate(encoding_type = fct_relevel(encoding_type,
                                         "flynet.only",
                                         "flynet.onoff",
                                         "flynet.minus.onoff",
                                         "alexnet.only",
                                         "alexnet.onoff",
                                         "alexnet.minus.onoff",
                                         "flynet.alexnet",
                                         "onoff.only"),
             encoding_family = case_match(encoding_type,
                                          c("flynet.only",
                                            "flynet.onoff",
                                            "flynet.minus.onoff") ~ "flynet",
                                          c("alexnet.only",
                                            "alexnet.onoff",
                                            "alexnet.minus.onoff") ~ "alexnet",
                                          c("flynet.alexnet") ~ "combined",
                                          c("onoff.only") ~ "onoff"),
             encoding_family = fct_relevel(encoding_family,
                                           "flynet",
                                           "alexnet",
                                           "combined"))
  }
  
  out <- perf %>% 
    ggplot(aes(x = encoding_type, y = perf)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = fold_num), alpha = 0.2) +
    geom_line(aes(group = 1), stat = "summary", fun = "mean") +
    geom_point(aes(color = encoding_family), size = 3, alpha = 0.2) +
    geom_pointrange(aes(color = encoding_family), size = 0.75, stat = "summary", fun.data = \(x) mean_se(x, mult = 1)) +
    guides(x = guide_axis(angle = 30)) +
    labs(x = "which encoding model?",
         y = sprintf("Predicted-actual BOLD correlation"),
         color = "Encoding model family",
         subtitle = glue::glue("Leave-one-subject-out cross-validation, N = {n_subjs}"))
  
  if (length(rois) > 1) {
    out <- out + facet_wrap(~ roi, scales = "free_y")
  }
  
  return (out)
}

# There isn't a pcor.combined target, but this expects two pcor targets 
# partial-ing the same two predictors, from amygdala and sc row-bound together
plot_encoding_performance_pcor <- function (perf_pcor_combined) {
  # these have subj_num because downstream of bold targets that have subject labels on the sub-mat files
  n_subjs <- length(unique(perf_pcor_combined$subj_num))
  
  perf_pcor_combined %>% 
    # this function does not know which encoding model names are coming into the pcor
    # so you need to relevel the order before you call this function
    ggplot(aes(x = term, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = subj_num), alpha = 0.2) +
    geom_point(alpha = 0.5) +
    geom_pointrange(stat = "summary", fun.data = "mean_se", color = "tomato") +
    facet_wrap(~ roi, scales = "free_y") +
    labs(x = "which encoding model?",
         y = "Predicted-observed BOLD PARTIAL correlation",
         subtitle = glue::glue("Leave-one-subject-out cross-validation, N = {n_subjs}"))
}

plot_encoding_selfreport_pcor <- function (summary_pcor) {
  summary_pcor %>% 
    mutate(vjust_label = case_match(term, 
                                    "scale(flynet)" ~ 0, 
                                    "scale(alexnet)" ~ 1),
           term = fct_recode(term, 
                             "looming" = "scale(flynet)", 
                             "object" = "scale(alexnet)"),
           rating_type = fct_recode(rating_type, "valence" = "pleasantness"),
           rating_type = fct_relevel(rating_type, "valence")) %>% 
    ggplot(aes(x = pcor_mean, y = fct_rev(rating_type), color = term)) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_pointrange(aes(xmin = pcor_ci95_lower, xmax = pcor_ci95_upper), 
                    position = position_dodge(width = 0.1)) + 
    geom_text(aes(label = term, vjust = vjust_label), 
              position = position_dodge(width = rel(0.4)),
              size = rel(5)) + 
    guides(color = "none") + 
    labs(x = sprintf("Partial \U03C1 \U00B1 95%% CI"), 
         y = NULL)
}

plot_encoding_controlled_flynet <- function (pattern_expressions) {
  pattern_expressions %>% 
    ggplot(aes(x = animal_type, y = pattern_expression, fill = direction)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0, linetype = "dotted")
}
