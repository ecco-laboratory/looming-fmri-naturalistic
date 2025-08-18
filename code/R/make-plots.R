## targets-optimized functions for generating plots from data targets ----

relabel_cols_for_plot_naturalistic <- function (data) {
  data %>% 
    mutate(loom_col = if_else(has_loom == 1, "Looming", "No looming"),
           animal_type = fct_relevel(animal_type, "food", "dog", "cat"))
}

relabel_cols_for_plot_controlled <- function (data) {
  data %>% 
    rename(loom_col = direction)
}

pivot_selfreport_longer <- function (data) {
  data %>% 
    pivot_longer(cols = starts_with("rating"), names_to = "construct", values_to = "rating", names_prefix = "rating_") %>% 
    mutate(construct = if_else(construct == "pleasantness", "valence (- to +)", construct),
           construct = fct_relevel(construct, "valence (- to +)"))
}

plot_selfreport_ratings <- function (beh_data) {
  n_subjs <- length(unique(beh_data$subj_num))
  
  preplot <- beh_data %>% 
    pivot_selfreport_longer()
  
  preplot %>% 
    # so that stuff will be plotted by subject
    group_by(construct, subj_num, loom_col, animal_type) %>% 
    summarize(rating = mean(rating), .groups = "drop") %>% 
    ggplot(aes(x = animal_type, y = rating, color = loom_col)) +
    # individual subjects spaghetti
    geom_line(aes(group = interaction(subj_num, loom_col)), alpha = 0.1) +
    # summary lines
    geom_line(aes(group = loom_col), stat = "summary", fun = "mean") +
    geom_pointrange(stat = "summary", 
                    fun.data = "mean_se") +
    labs(x = "Object type",
         y = "Self-report rating",
         color = "Motion",
         subtitle = glue::glue("N = {n_subjs} participants")) +
    facet_wrap(~ construct)
  
}

plot_selfreport_ratings_by_stim <- function (beh_data) {
  beh_data %>% 
    pivot_selfreport_longer() %>% 
    group_by(construct, video_id, loom_col, animal_type) %>% 
    summarize(rating_mean = mean(rating), 
              rating_se = sd(rating)/sqrt(n())) %>% 
    group_by(construct) %>% 
    # so that x is facet-specific
    mutate(rating_rank = row_number(rating_mean),
           rating_rank = rating_rank/max(rating_rank)) %>% 
    ggplot(aes(x = rating_rank, y = rating_mean, color = animal_type, shape = loom_col)) + 
    geom_pointrange(aes(ymin = rating_mean - rating_se, ymax = rating_mean + rating_se), fatten = 2) + 
    scale_shape_manual(values = c(19, 3)) + 
    facet_wrap(facets = vars(construct), ncol = 1) +
    labs(x = NULL,
         y = "Mean self-report rating across subjects ± 1 SE",
         color = "Object type",
         shape = "Motion")
}

plot_norm_ratings <- function (norm_data) {
  preplot <- norm_data %>% 
    pivot_selfreport_longer()
  
  preplot %>% 
    ggplot(aes(x = animal_type, y = rating, color = loom_col)) +
    # points, not individual subjects spaghetti, because we don't know that every norming participant saw every stimulus
    geom_jitter(alpha = 0.1, width = 0.1) +
    # summary lines
    geom_line(aes(group = loom_col), stat = "summary", fun = "mean") +
    geom_pointrange(stat = "summary", 
                    fun.data = "mean_se") +
    labs(x = "Object type",
         y = "Self-report rating",
         color = "Motion") +
    facet_wrap(~ construct)
}

# will use label_both labeller to put the name of the facet col before the unit/voxel num
plot_sample_timecourse <- function (plot_data, y_col, facet_col, facet_title_pos = "left", line_color = "black") {

  plot_out <- ggplot(plot_data,
         aes(x = tr_num, y = {{y_col}})) + 
    geom_line(color = line_color, linewidth = 0.5) + 
    labs(x = "TR (concatenated across runs)")
  
  # only facet if there would be more than one panel
  if (length(unique(pull(plot_data, {{facet_col}}))) > 1) {
    plot_out <- plot_out +
      facet_wrap(facets = vars({{facet_col}}),
                 ncol = 1,
                 scales = "free_y",
                 strip.position = facet_title_pos,
                 labeller = \(x) label_both(x, sep = " "))
  }
  
  return (plot_out)
}

# takes a named list of targets of cross-validated encoding predictions, a la pred.encoding.xval_{ROI}
# the list names will become the facet labels!
plot_sample_timecourse_encoding <- function (list_paths_pred_encoding, subj = 1) {
  # only 1 subject at a time!
  stopifnot(length(subj) == 1)
  
  encoding_names <- names(list_paths_pred_encoding)
  
  preplot <- list_paths_pred_encoding %>% 
    map(\(x) read_csv(x, col_names = FALSE)) %>% 
    bind_rows(.id = "encoding_type") %>% 
    select(encoding_type, subj_num = X1, everything()) %>% 
    # just do it for subject 1 by default because this is for a schematic
    # but set it as an arg so you COULD set it later if you needed to
    filter(subj_num == subj) %>% 
    # to relabel and reorder them into the names and level order specified by the argument. for plot facet order
    mutate(encoding_type = fct_relevel(encoding_type, !!encoding_names))
  
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

# to plot selected units from a single subject/single encoding model target
plot_sample_timecourse_activation <- function (target_subj_activation, units_to_plot, line_color) {
  target_subj_activation %>% 
    map(\(x) read_csv(x, col_names = FALSE)) %>% 
    bind_rows() %>% 
    mutate(tr_num = 1:n()) %>% 
    # that's also why the pivot is inside this map
    pivot_longer(cols = -tr_num, 
                 names_to = "unit", 
                 values_to = "activation",
                 names_transform = list(unit = \(x) as.integer(str_sub(x, start = 2L)))) %>% 
    filter(unit %in% units_to_plot) %>%
    plot_sample_timecourse(y_col = activation, facet_col = unit, line_color = line_color) +
    labs(y = "neural network model activation")
}

# takes a named list of 2 across-run single-subject model activation targets
# bc the other schematics also facet by encoding model type
plot_sample_timecourse_activation_combined <- function (targets_subj_activation) {
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
plot_sample_timecourse_bold <- function (path_bold_masked_mat, voxels_to_plot, line_color = "black") {
  path_bold_masked_mat %>% 
    R.matlab::readMat() %>% 
    pluck("DATA") %>% 
    as_tibble(.name_repair = "unique") %>% 
    mutate(tr_num = 1:n()) %>% 
    pivot_longer(cols = starts_with("..."), 
                 names_to = "voxel", 
                 values_to = "bold", 
                 names_transform = list(voxel = \(x) as.integer(str_sub(x, start = 4L)))) %>% 
    # again to cut down the number of voxels getting plotted
    filter(voxel %in% voxels_to_plot) %>% 
    plot_sample_timecourse(y_col = bold, facet_col = voxel, facet_title_pos = "bottom", line_color = line_color) +
    labs(y = "BOLD")
}

# the betas targets are allsubs, so this filters for one
plot_sample_encoding_betas <- function (path_betas, subj = 1, viridis_option = "D") {
    read_csv(path_betas, col_names = FALSE) %>% 
    filter(X1 == subj) %>% 
    # the first row is intercept
    slice(-1) %>% 
    select(-X1) %>% 
    mutate(unit_num = 1:n()) %>% 
    pivot_longer(cols = -unit_num, 
                 names_to = "voxel_num", 
                 values_to = "beta", 
                 names_transform = list(voxel_num = \(x) as.integer(str_sub(x, start = 2L)) - 1)) %>% 
    ggplot(aes(x = voxel_num, y = unit_num, fill = beta)) + 
    geom_raster() + 
    # so that more units will make a taller plot
    coord_fixed(ratio = 1/4) +
    guides(fill = "none") +
    scale_fill_viridis_c(option = viridis_option) +
    labs(x = "voxel", y = "neural network unit")
}

plot_confusion <- function(object_preds_unnested, 
                           brewer_option,
                           label_family = "Nimbus Sans", 
                           facet_var_row = NULL, 
                           facet_var_col = NULL, 
                           facet_var_wrap = NULL) {
  facet_vars <- enquos(facet_var_row, facet_var_col, facet_var_wrap)
  facet_var_row <- enquo(facet_var_row)
  facet_var_col <- enquo(facet_var_col)
  facet_var_wrap <- enquo(facet_var_wrap)
  
  plot_data <- object_preds_unnested %>% 
    count(!!!facet_vars, .obs, .pred) %>% 
    # bc complete uses all implicit factor levels, which may re-input intentionally empty categories like food
    mutate(across(c(.obs, .pred), fct_drop)) %>%
    complete(!!!facet_vars, .obs, .pred, fill = list(n = 0)) %>%
    group_by(pick(!!!facet_vars, .obs)) %>% 
    mutate(prob = n/sum(n)) %>% 
    ungroup() %>% 
    mutate(diagonal_outline_width = if_else(.pred == .obs, .5, NA))
  
  plot_out <- plot_data %>% 
    ggplot(aes(x = .pred, y = fct_rev(.obs))) + 
    geom_raster(aes(fill = prob)) +
    # 2025-07-22: this line would draw an outline the diagonal cells 
    # but I think we don't need it now that on-diagonal performance is pretty good
    # geom_tile(aes(linewidth = diagonal_outline_width), width = 1, height = 1, alpha = 0, color = "grey85") + 
    geom_text(aes(label = round(prob, 2)), color = "black", family = label_family) + 
    scale_fill_distiller(palette = brewer_option, direction = 1) +
    # 2025-07-22: stop putting a fill legend to save space bc prob is already written right onto the cells
    guides(fill = "none") +
    labs(x = "Predicted class", y = "True class")
  
  if (!quo_is_null(facet_var_row) | !quo_is_null(facet_var_col)) {
    plot_out <- plot_out +
      facet_grid(rows = vars(!!facet_var_row), cols = vars(!!facet_var_col))
  } else if (!quo_is_null(facet_var_wrap)) {
    plot_out <- plot_out +
      facet_wrap(facets = vars(!!facet_var_wrap))
  }

  
  return (plot_out)
}

plot_auroc <- function (object_preds, grouping_col = NULL, include_overall = TRUE) {
  grouping_col <- enquo(grouping_col)
  
  preplot_data <- object_preds %>% 
    select(-coefs) %>% 
    unnest(preds) %>% 
    # 2025-07-14 at this point food is being excluded from all the main analyses
    filter(animal_type != "food")
  
  preplot_data_wider <- preplot_data %>%
    unnest_wider(.preds) %>% 
    # if it's there we don't want it srsly!!!
    select(-any_of(".pred_outcome_food")) %>% 
    mutate(.obs = fct_drop(.obs))
  
  plot_data <- preplot_data_wider %>%
    roc_curve(truth = .obs, starts_with(".pred_outcome"))
  
  if (!quo_is_null(grouping_col)) {
    grouping_col_name <- as_name(grouping_col)
    
    if (include_overall) {
      plot_data_grouped <- preplot_data_wider %>%
        group_by(pick({{grouping_col}})) %>% 
        roc_curve(truth = .obs, starts_with(".pred_outcome"))
      
      plot_data %<>%
        mutate(!!grouping_col_name := "overall") %>% 
        bind_rows(plot_data_grouped) %>% 
        mutate(!!grouping_col_name := fct_relevel(!!grouping_col, "overall", after = Inf))
      
      if (grouping_col_name == "animal_type") {
        plot_data %<>%
          mutate(!!grouping_col_name := fct_relevel(!!grouping_col, "dog", "cat", "frog", "spider"))
      }
    } else {
      plot_data <- plot_data_grouped
    }
    
    # not using autoplot because if you tag the overall back onto the grouped
    # it loses the roc_df class which autoplot needs to choose the correct method
    out <- plot_data %>% 
      ggplot(aes(x = 1 - specificity, y = sensitivity, color = animal_type)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") + 
      geom_path()
    
    return (out)
  } else if (".level" %in% names(plot_data)) {
    # if it's from a multiclass roc_curve, the class level comes in direct from roc_curve
    # and need to calculate the overall separately, bc the threshold values aren't consistent so can't group_by summarize
    if (include_overall) {
      plot_data_overall <- preplot_data %>% 
        unnest_longer(.preds, values_to = ".pred_value", indices_to = "animal_to_guess") %>% 
        mutate(animal_to_guess = str_remove(animal_to_guess, ".pred_outcome_"), 
               .obs_binary = if_else(animal_type == animal_to_guess, "yes", "no"), 
               .obs_binary = factor(.obs_binary, levels = c("yes", "no"))) %>% 
        roc_curve(truth = .obs_binary, .pred_value) %>% 
        mutate(.level = "overall")
      
      plot_data %<>%
        bind_rows(plot_data_overall)
    }
    
    out <- plot_data %>% 
      # will put overall last if there. can do it now that the roc curves have been calculated
      mutate(.level = fct_relevel(.level, "dog", "cat", "frog", "spider")) %>% 
      ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") + 
      geom_path()
    
    return (out)
  } else {
    return (autoplot(plot_data))
  }
  
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


# since we now have to do enough of these that I may as well standardize the appearance
plot_spaghetti_by_model <- function (data, model_col, y_col, subj_col, color_col = NULL, spaghetti_alpha = 0.2) {
  # must enquo no matter what so that you can test the type without evaluating the var name
  color_col <- enquo(color_col)
  if (quo_is_null(color_col)) color_col <- enquo(model_col)
  
 data %>% 
    ggplot(aes(x = {{model_col}}, y = {{y_col}})) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = {{subj_col}}), alpha = spaghetti_alpha) +
    geom_line(aes(group = 1), stat = "summary", fun = "mean") +
    geom_point(aes(color = {{color_col}}), size = 3, alpha = spaghetti_alpha) +
    geom_pointrange(aes(color = {{color_col}}), size = 0.75, stat = "summary", fun.data = \(x) mean_se(x, mult = 1))
}

plot_encoding_performance <- function (perf_combined, encoding_types = NULL) {
  n_subjs <- length(unique(perf_combined$fold_num))
  
  plot_data <- perf_combined
  
  if (!is.null(encoding_types)) {
    plot_data %<>%
      filter(encoding_type %in% encoding_types) %>% 
      mutate(encoding_type = fct_relevel(encoding_type, !!encoding_types),
             encoding_type = fct_recode(encoding_type, !!!encoding_types),
             encoding_family = encoding_type)
  } else {
    plot_data %<>%
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
  
  out <- plot_data %>% 
    plot_spaghetti_by_model(model_col = encoding_type,
                            y_col = perf,
                            subj_col = fold_num,
                            color_col = encoding_family) +
    guides(x = guide_axis(angle = 30)) +
    labs(x = "which encoding model?",
         y = "Performance (Pearson's _r_)",
         color = "Encoding model family",
         subtitle = glue::glue("Leave-one-subject-out cross-validation, N = {n_subjs}"))
  
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

# this one works on the regular correlation summary OR the pcor summary! but not both together, that requires more processing
plot_encoding_selfreport <- function (summary_encoding_selfreport, cor_col_prefix = "pcor", xlab = "Partial ρ ± 95% CI") {
  summary_encoding_selfreport %>% 
    rename_with(\(x) str_replace(x, cor_col_prefix, "cor"), starts_with(cor_col_prefix)) %>% 
    mutate(vjust_label = case_match(model_type, 
                                    "flynet" ~ 0, 
                                    "alexnet" ~ 1),
           model_type = fct_recode(model_type, 
                             "looming" = "flynet", 
                             "object" = "alexnet"),
           rating_type = fct_recode(rating_type, "valence" = "pleasantness"),
           rating_type = fct_relevel(rating_type, "valence")) %>% 
    ggplot(aes(x = cor_mean, y = fct_rev(rating_type), color = model_type)) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_pointrange(aes(xmin = cor_ci95.lower, xmax = cor_ci95.upper), 
                    position = position_dodge(width = 0.1)) + 
    geom_text(aes(label = model_type, vjust = vjust_label), 
              position = position_dodge(width = rel(0.4)),
              size = rel(5)) + 
    guides(color = "none") + 
    labs(x = xlab, 
         y = NULL)
}

plot_encoding_selfreport_combined <- function (summary_encoding_selfreport_cor,
                                               summary_encoding_selfreport_pcor) {
  # rcor for "regular correlation". dumb but I just don't want to risk overloading cor() in the namespace
  # harmonize the colnames before row-binding
  bind_rows("full" = summary_encoding_selfreport_cor %>% 
              rename_with(\(x) str_replace(x, "correlation", "cor"), everything()),
            "partial" = summary_encoding_selfreport_pcor %>% 
              rename_with(\(x) str_replace(x, "pcor", "cor"), everything()),
            .id = "cor_type") %>% 
    mutate(model_type = fct_recode(model_type, 
                                   "looming" = "flynet", 
                                   "object" = "alexnet"),
           rating_type = fct_recode(rating_type, "valence" = "pleasantness"),
           rating_type = fct_relevel(rating_type, "valence")) %>% 
    ggplot(aes(x = cor_type, y = cor_mean, color = model_type, fill = model_type)) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    geom_ribbon(aes(ymin = cor_ci95.lower, ymax = cor_ci95.upper, group = model_type), alpha = 0.4, linewidth = 0) +
    geom_line(aes(group = model_type)) +
    geom_point(size = 3) + 
    labs(x = "correlation type", 
         y = "correlation",
         color = "Feature type",
         fill = "Feature type") +
    facet_wrap(~ rating_type)
}

plot_selfreport_statistics <- function (values,
                                        x_col,
                                        y_col,
                                        color_col = NULL,
                                        color_constant = "black") {
  
  plot_data <- values %>% 
    mutate(rating_type = if_else(rating_type == "pleasantness", "valence (- to +)", rating_type),
           rating_type = fct_relevel(rating_type, "valence (- to +)"))
  
  if (!quo_is_null(enquo(color_col))) {
    out <- plot_data %>% 
      ggplot(aes(x = {{x_col}}, y = {{y_col}}, color = {{color_col}})) + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      geom_jitter(alpha = 0.1, width = 0.1) + 
      geom_pointrange(stat = "summary", fun.data = "mean_se")
  } else {
    out <- plot_data %>% 
      ggplot(aes(x = {{x_col}}, y = {{y_col}})) + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      geom_jitter(alpha = 0.1, width = 0.1, color = color_constant) + 
      geom_pointrange(stat = "summary", fun.data = "mean_se", color = color_constant)
  }
  
  return (out)
  
}

plot_encoding_controlled_flynet <- function (pattern_expressions, y_var, y_label) {
  pattern_expressions %>% 
    ggplot(aes(x = animal_type, y = {{y_var}}, color = direction)) + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(aes(group = interaction(subj_num, direction)), alpha = 0.2) + 
    geom_line(aes(group = direction), stat = "summary", fun = "mean") + 
    geom_pointrange(stat = "summary", fun.data = "mean_se") +
    labs(x = "Object type",
         y = y_label,
         color = "Controlled motion\ndirection")
}

plot_cormats_encoding_object_selfreport <- function (beta.comparison) {
  
  beta.comparison %>% 
    mutate(across(c(row, col), 
                  \(x) if_else(x == "outcome_loom",
                               "looming_looming",
                               x) %>% 
                    str_replace("outcome_", "object_object: ") %>% 
                    str_remove("rating_") %>% 
                    # unicode right arrow, for looming -> self-report: etc
                    str_replace("_", "→") %>% 
                    fct_relevel("looming→looming", 
                                "object→object: dog", 
                                "object→object: cat", 
                                "object→object: frog", 
                                "object→object: spider",
                                "looming→unpleasantness",
                                "looming→arousal",
                                "looming→fear",
                                "object→unpleasantness",
                                "object→arousal",
                                "object→fear")),
           ci_excludes_0 = sign(cor_ci95.lower) == sign(cor_ci95.upper)) %>% 
    ggplot(aes(x = col, y = fct_rev(row), fill = cor_true)) + 
    geom_raster() + 
    geom_text(aes(label = round(cor_true, 2),
                  color = ci_excludes_0)) + 
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey64")) +
    scale_fill_gradient2() +
    guides(x = guide_axis(angle = 30),
           color = "none") +
    labs(x = NULL, y = NULL, fill = "Correlation")
}

plot_parcel_connectivity_ggseg <- function (parcels_long,
                                            fill_col,
                                            pval_col,
                                            ggseg_sides,
                                            ggseg_position = "horizontal",
                                            positive_only = TRUE,
                                            max_fill_tval = 5,
                                            threshold_p_outline = .05,
                                            viridis_palette = "viridis") {
  if (positive_only) {
    plot_data <- parcels_long %>% 
      mutate(this_fill_col = pmax({{fill_col}}, 0),
             this_outline_col = {{pval_col}} < threshold_p_outline & this_fill_col > 0)
    
    this_scale_fill <- scale_fill_viridis_c(limits = c(0, max_fill_tval), oob=scales::squish, option = viridis_palette)
  } else {
    plot_data <- parcels_long %>% 
      rename(this_fill_col = {{fill_col}}) %>% 
      mutate(this_outline_col = {{pval_col}} < threshold_p_outline)
    this_scale_fill <- scale_fill_gradient2(limits = c(-max_fill_tval, max_fill_tval), oob=scales::squish)
  }
  
  plot_data %>% 
    relabel_glasser_clt2ggseg() %>% 
    brain_join(ggsegGlasser::glasser) %>% 
    filter(side %in% ggseg_sides) %>% 
    reposition_brain(ggseg_position) %>% 
    arrange(this_fill_col) %>%
    ggplot() +
    geom_sf(aes(fill = this_fill_col, color = this_outline_col), linewidth = 0.5) +
    geom_sf_label(aes(label = if_else(this_outline_col, region, NA_character_)), size = 3) +
    this_scale_fill +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), na.translate = FALSE) +
    guides(color = "none") +
    labs(fill = "t-value")
}

plot_parcel_connectivity_ggseg_cat <- function (parcels_long,
                                                fill_col,
                                                outline_col,
                                                ggseg_sides,
                                                ggseg_position = "horizontal",
                                                threshold_p_outline = .05) {
  parcels_long %>% 
    relabel_glasser_clt2ggseg() %>% 
    brain_join(ggsegGlasser::glasser) %>% 
    filter(side %in% ggseg_sides) %>% 
    reposition_brain(ggseg_position) %>% 
    arrange(desc({{outline_col}})) %>%
    ggplot() +
    geom_sf(aes(fill = {{fill_col}}, color = {{outline_col}}), linewidth = 0.5) +
    geom_sf_label(aes(label = if_else({{outline_col}} == "apriori", region, NA_character_)), size = 3) +
    guides(color = "none") +
    labs(fill = "which is significant?")
}

plot_parcel_connectivity_scatter <- function (parcel_values_relabeled,
                                              facet_col) {
  facet_col <- enquo(facet_col)
  
  plot_out <- parcel_values_relabeled %>% 
    plot_spaghetti_by_model(model_col = model_type, y_col = value, subj_col = fold_num, spaghetti_alpha = 0.1) +
    facet_wrap(facets = vars(!!facet_col), nrow = 2, strip.position = "left") +
    guides(color = "none", x = guide_axis(angle = 30)) +
    labs(x = "Connectivity type",
         y = "Mean connectivity (Pearson's _r_)")

  return (plot_out)
}

plot_sig_sim <- function (sig_sim) {
  plot_data <- sig_sim %>% 
    pivot_longer(cols = -c(model_type, fold_num), names_to = "signature", values_to = "similarity") %>% 
    mutate(signature = fct_rev(fct_reorder(signature, similarity, .fun = mean, .na_rm = TRUE)))
  
  plot_data %>% 
    ggplot(aes(x = signature, y = similarity)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(aes(color = model_type), 
               position = position_dodge(width = 0.5), 
               size = 3, 
               alpha = 0.1) +
    geom_pointrange(aes(color = model_type), 
                    position = position_dodge(width = 0.5), 
                    size = 0.75, 
                    stat = "summary", fun.data = "mean_se") +
    guides(x = guide_axis(angle = 30)) +
    labs(x = "Whole-brain signature",
         y = "Pattern expression (cosine similarity)",
         color = "Connectivity type")
}
