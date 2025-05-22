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
         color = "Motion",
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

plot_encoding_object_acc <- function (acc_summary, y_description, encoding_labels = NULL) {
  
  if (!is.null(encoding_labels)) {
    acc_summary %<>%
      mutate(encoding_type = fct_relevel(encoding_type, !!encoding_labels),
             encoding_type = fct_recode(encoding_type, !!!encoding_labels))
  }
  
  acc_summary %>% 
    filter(.metric == "accuracy") %>% 
    mutate(pval_text = if_else(pval < .001, "p < .001", glue::glue("p = {round(pval, 3)}"))) %>% 
    ggplot(aes(x = encoding_type, y = estimate, fill = encoding_type)) + 
    geom_col() +
    geom_text(aes(label = pval_text), vjust = 0, nudge_y = .01) +
    guides(fill = "none") +
    labs(x = "which encoding model?",
         y = glue::glue("Classification accuracy\n({y_description})"))
}

plot_encoding_object_acc_by_object <- function (acc_summary, y_description, encoding_labels = NULL) {
  
  if (!is.null(encoding_labels)) {
    acc_summary %<>%
      mutate(encoding_type = fct_relevel(encoding_type, !!encoding_labels),
             encoding_type = fct_recode(encoding_type, !!!encoding_labels))
  }
  
  acc_summary %>% 
    filter(.metric == "accuracy") %>% 
    mutate(animal_type = fct_relevel(animal_type, "dog", "cat", "frog", "spider"),
           pval_text = if_else(pval < .001, "p < .001", glue::glue("p = {round(pval, 3)}"))) %>% 
    ggplot(aes(x = animal_type, y = estimate, fill = encoding_type)) + 
    geom_col() +
    geom_text(aes(label = pval_text), vjust = 0, nudge_y = .01) +
    guides(fill = "none") +
    labs(x = "Object category",
         y = glue::glue("Classification accuracy\n({y_description})")) +
    facet_wrap(~ encoding_type)
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
    geom_pointrange(aes(xmin = pcor_ci95.lower, xmax = pcor_ci95.upper), 
                    position = position_dodge(width = 0.1)) + 
    geom_text(aes(label = term, vjust = vjust_label), 
              position = position_dodge(width = rel(0.4)),
              size = rel(5)) + 
    guides(color = "none") + 
    labs(x = sprintf("Partial \U03C1 \U00B1 95%% CI"), 
         y = NULL)
}

plot_predicted_alexnet_categories <- function (path_alexnet_activations, stim_labels, path_imagenet_categories, lump_prop = .1) {
  read_csv(path_alexnet_activations) %>% 
    right_join(stim_labels, 
               by = "video") %>% 
    select(-frame) %>% 
    nest(.by = c(video, animal_type, looming), .key = "activations") %>% 
    mutate(activations = map(activations, as.matrix)) %>% 
    mutate(max_indices_from0 = map(activations, \(x) apply(x, 1, which.max) - 1)) %>% 
    select(-activations) %>% 
    unchop(max_indices_from0) %>% 
    left_join(read_csv(path_imagenet_categories), by = c("max_indices_from0" = "index")) %>% 
    group_by(animal_type, looming) %>% 
    mutate(categories_lumped = fct_lump_prop(categories, prop = lump_prop, other_level = "other")) %>% 
    count(animal_type, looming, categories_lumped) %>% 
    arrange(animal_type, looming, desc(categories_lumped)) %>% 
    mutate(prop = n/sum(n),
           bar_max = cumsum(prop), 
           bar_min = coalesce(lag(bar_max), 0), 
           bar_mid = bar_min + (bar_max-bar_min)/2) %>% 
    ggplot(aes(x = looming, y = prop, fill = categories_lumped)) + 
    geom_col() + 
    geom_text(aes(y = bar_mid, label = categories_lumped)) + 
    facet_wrap(~animal_type) + 
    scale_fill_brewer(type = "qual", palette = "Set3") +
    guides(fill = "none") +
    labs(x = NULL,
         y = "Top-1 ImageNet superordinate class probability")
}

plot_encoding_controlled_flynet <- function (pattern_expressions) {
  pattern_expressions %>% 
    ggplot(aes(x = animal_type, y = pattern_expression, color = direction)) + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = interaction(subj_num, direction)), alpha = 0.2) + 
    geom_line(aes(group = direction), stat = "summary", fun = "mean") + 
    geom_pointrange(stat = "summary", fun.data = "mean_se") +
    labs(x = "Object type",
         y = "Looming pattern expression strength (AU)",
         color = "Controlled motion\ndirection")
}

plot_cormats_encoding_object_selfreport <- function (beta.comparison_flynet,
                                                     beta.comparison_alexnet,
                                                     remove_food = TRUE) {
  
  plot_data <- bind_rows("Looming encoding model" = beta.comparison_flynet,
            "Object encoding model" = beta.comparison_alexnet,
            .id = "encoding_type")
  
  if (remove_food) {
    plot_data %<>%
      filter(row != "outcome_food", col != "outcome_food")
  }
  
  plot_data %>% 
    mutate(across(c(row, col), 
                  \(x) if_else(x == "outcome_loom",
                               "looming",
                               x) %>% 
                    str_replace("outcome_", "object: ") %>% 
                    str_replace("rating_", "self-report: ") %>% 
                    fct_relevel("looming", 
                                "object: dog", 
                                "object: cat", 
                                "object: frog", 
                                "object: spider", 
                                "self-report: unpleasantness")),
           ci_excludes_0 = sign(cor_ci95.lower) == sign(cor_ci95.upper)) %>% 
    ggplot(aes(x = col, y = fct_rev(row), fill = cor_true)) + 
    geom_raster() + 
    geom_text(aes(label = round(cor_true, 2),
                  color = ci_excludes_0)) + 
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey64")) +
    scale_fill_gradient2() +
    guides(x = guide_axis(angle = 30),
           color = "none") +
    facet_wrap(~ encoding_type) +
    labs(x = NULL, y = NULL, fill = "Correlation")
}
