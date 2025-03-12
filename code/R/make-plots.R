## targets-optimized functions for generating plots from data targets ----

plot_selfreport_ratings <- function (beh_data) {
  n_subjs <- length(unique(beh_data$subj_num))
  
  preplot <- beh_data %>% 
    mutate(has_loom = if_else(has_loom == 1, "Looming", "No looming"),
           animal_type = fct_relevel(animal_type, "food", "dog", "cat")) %>% 
    pivot_longer(cols = ends_with("rating"), names_to = "construct", values_to = "rating") %>% 
    mutate(construct = str_remove(construct, "_rating"),
           construct = if_else(construct == "pleasantness", "valence (- to +)", construct),
           construct = fct_relevel(construct, "valence (- to +)"))
  
  preplot %>% 
    # so that stuff will be plotted by subject
    group_by(subj_num, animal_type, has_loom, construct) %>%
    summarize(rating = mean(rating), .groups = "drop") %>%
    ggplot(aes(x = animal_type, y = rating, color = has_loom)) +
    # individual subjects spaghetti
    geom_line(aes(group = interaction(subj_num, has_loom)), alpha = 0.2) +
    # summary lines
    geom_line(aes(group = has_loom), stat = "summary", fun = "mean") +
    geom_pointrange(stat = "summary", 
                    fun.data = "mean_se") +
    labs(x = "Object type",
         y = "Self-report rating",
         color = NULL,
         subtitle = glue::glue("N = {n_subjs} participants")) +
    facet_wrap(~ construct) +
    this_theme
  
}

# takes a long df target of cross-validated encoding predictions, a la pred.encoding.xval_{ROI}
plot_sample_timecourse_encoding <- function (pred_encoding_xval, encoding_models = NULL, fold = 1) {
  # only 1 subject at a time!
  stopifnot(length(fold) == 1)
  
  preplot <- pred_encoding_xval %>% 
    rename(fold_num = X1) %>% 
    select(encoding_type, fold_num, everything()) %>% 
    # just do it for subject 1 by default because this is for a schematic
    # but set it as an arg so you COULD set it later if you needed to
    filter(fold_num == fold)
  
  if (!is.null(encoding_models)) {
    preplot %<>%
      filter(encoding_type %in% encoding_models) %>% 
      # to reorder them into the level specified by the argument. for plot facet order
      mutate(encoding_type = factor(encoding_type, levels = encoding_models))
  }
  
  preplot %<>%
    group_by(encoding_type) %>% 
    mutate(tr_num = 1:n()) %>% 
    ungroup() %>% 
    pivot_longer(cols = starts_with("X"), 
                 names_to = "voxel_num", 
                 values_to = "bold_pred", 
                 names_transform = list(voxel_num = \(x) as.integer(str_sub(x, start = 2L)) - 1L)) %>%
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

plot_encoding_performance <- function (perf_combined) {
  n_subjs <- length(unique(perf_combined$fold_num))
  
  perf_combined %>% 
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
                                         "combined")) %>% 
    ggplot(aes(x = encoding_type, y = perf)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = fold_num), alpha = 0.2) +
    geom_point(alpha = 0.5) +
    geom_pointrange(aes(color = encoding_family), stat = "summary", fun.data = "mean_se") +
    guides(x = guide_axis(angle = 30)) +
    facet_wrap(~ roi, scales = "free_y") +
    labs(x = "which encoding model?",
         y = "Predicted-observed BOLD timecourse correlation",
         color = "Encoding model family",
         subtitle = glue::glue("Leave-one-subject-out cross-validation, N = {n_subjs}"))
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
