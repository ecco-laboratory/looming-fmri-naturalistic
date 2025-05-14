make_targets_fmri_by.run <- function (n_runs, task = "controlled", additional_targets = NULL) {
  
  if (task == "controlled") {
    subtarget_events <- tar_target(name = events,
                                   # n_trs fed in here is not including any TRs during the 8-second wait period
                                   # it counts from the first TR whose acquisition crosses the 8-second boundary
                                   command = parse_events_nback(events.raw, 
                                                                tr_duration = task_defaults_list$tr_duration, 
                                                                n_trs = task_defaults_list$n_trs_kept))
    
  } else if (task == "naturalistic") {
    subtarget_events <- tar_target(name = events,
                                   # n_trs fed in here is not including any TRs during the 8-second wait period
                                   # it counts from the first TR whose acquisition crosses the 8-second boundary
                                   command = parse_events_naturalistic(events.raw, 
                                                                       conditions_base = stims,
                                                                       tr_duration = task_defaults_list$tr_duration, 
                                                                       n_trs = task_defaults_list$n_trs_kept))
  }
  tar_map(
    values = tibble(run_num = 1:n_runs) %>% 
      mutate(run_bids = sprintf("run-%02d", run_num),
             task = task,
             task_bids = paste0("task-", task)),
    tar_target(name = events.raw,
               command = get_raw_events(subject, task_bids, run_bids),
               format = "file"),
    subtarget_events,
    tar_target(name = confounds,
               command = get_raw_confounds(subject, task_bids, run_bids),
               format = "file"),
    tar_target(name = confounds.prespm,
               command = get_confounds_motion(file_path = confounds,
                                              # fyi again these are the same for both tasks
                                              tr_duration = task_defaults_list$tr_duration,
                                              disdaq_duration = task_defaults_list$disdaq_duration),
               format = "file"),
    # defining these separately instead of in a tar_eval together 
    # because we need them to be accessible as components of the list tar_map
    tar_target(name = events.prespm.boxcar,
               command = format_events_matlab(onsets = events,
                                              raw_path = events.raw,
                                              onset_type = "boxcar"),
               format = "file"),
    tar_target(name = events.prespm.endspike,
               command = format_events_matlab(onsets = events,
                                              raw_path = events.raw,
                                              onset_type = "endspike"),
               format = "file"),
    tar_target(name = bold_gz,
               command = get_bold_gz(subject, task_bids, run_bids),
               format = "file"),
    tar_target(name = bold,
               command = gunzip_bold(bold_gz),
               format = "file"),
    tar_target(name = bold.smoothed,
               command = spm_smooth(bold_path = bold,
                                    # THIS TAKES THE FULL, UNDISCARDED NUMBER OF TRs
                                    # it seems more future/forget-proof not to shorten the smoothed timeseries relative to the base timeseries
                                    n_trs = task_defaults_list$n_trs,
                                    # RIGHT NOW WE ARE ONLY SMOOTHING AT A 4 MM KERNEL! LOOKEE HERE
                                    kernel = 4,
                                    script = matlab_spmbatch_smooth),
               format = "file"),
    additional_targets,
    names = run_bids
  )
}

make_targets_fmri_by.subject <- function(df_participants, targets_by.run, contrast_names, task = "controlled", additional_targets = NULL) {
  
  task_bids <- paste0("task-", task)
  
  tar_map(
    values = df_participants %>% 
      mutate(task = task,
             task_bids = paste0("task-", task)),
    targets_by.run,
    tar_combine(name = bold.smoothed,
                targets_by.run[["bold.smoothed"]]),
    tar_combine(name = events.prespm.boxcar,
                targets_by.run[["events.prespm.boxcar"]]),
    tar_combine(name = events.prespm.endspike,
                targets_by.run[["events.prespm.endspike"]]),
    tar_combine(name = confounds.prespm,
                targets_by.run[["confounds.prespm"]]),
    ### level 1 SPM.mat models ----
    tar_target(name = spm.level1.boxcar,
               # the output file to be tracked is the SPM.mat file
               command = spm_spec_est_level1(model_path = file.path("ignore", "models", task_bids, "acq-mb8", subject, "model-boxcar", "smoothed-4mm"),
                                             tr_duration = task_defaults_list$tr_duration,
                                             trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                             # these three need to take vectors containing the values for each run
                                             bolds = bold.smoothed,
                                             events_prespm = events.prespm.boxcar,
                                             confounds_prespm = confounds.prespm,
                                             script = matlab_spmbatch_spec_est_level1),
               format = "file"),
    tar_target(name = spm.level1.endspike,
               # the output file to be tracked is the SPM.mat file
               command = spm_spec_est_level1(model_path = file.path("ignore", "models", task_bids, "acq-mb8", subject, "model-endspike", "smoothed-4mm"),
                                             tr_duration = task_defaults_list$tr_duration,
                                             trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                             # these three need to take vectors containing the values for each run
                                             bolds = bold.smoothed,
                                             events_prespm = events.prespm.endspike,
                                             confounds_prespm = confounds.prespm,
                                             script = matlab_spmbatch_spec_est_level1),
               format = "file"),
    ### readably named contrast targets by subject ----
    # defining contrasts individually (CRY) so that we can keep them accessible through one level of tar_map nesting
    make_contrast_targets_by.subject(spm.level1.boxcar, contrast_names, "boxcar"),
    make_contrast_targets_by.subject(spm.level1.endspike, contrast_names, "endspike"),
    additional_targets,
    names = subject
  )
}

# use tar_target_raw to construct a list of readably named targets for each level 1 contrast
make_contrast_targets_by.subject <- function (level1, contrast_names, model_type) {
  # use ensym() instead of enquo() because the tar_target_raw defused command 
  # just wants symbol names that will get evaluated at some miraculous point
  level1 <- ensym(level1)

  out <- list()
  for (i in 1:length(contrast_names)) {
    out <- c(out,
             tar_target_raw(name = paste("con", contrast_names[i], model_type, sep = "."),
                            command = expr(get_spm_level1_contrast(!!level1, !!i)),
                            format = "file"))
  }
  return (out)
}

make_targets_fmri_across.subject <- function (targets_by.subject, contrast_names, task = "controlled", additional_targets = NULL) {
  
  list(
    targets_by.subject,
    tar_eval(
      expr = tar_combine_raw(name = target_name,
                             targets_by.subject[[target_name]]),
      values = crossing(target_name = contrast_names,
                        model_type = c("boxcar", "endspike")) %>% 
        mutate(target_name = paste("con", target_name, model_type, sep = "."))
    ),
    tar_eval(
      tar_target(
        name = target_name,
        command = spm_spec_est_level2(model_path = file.path("ignore", "models", task_bids, "acq-mb8", "group", paste0("model-", model_type), "smoothed-4mm", contrast),
                                      # this need to take a vector containing the values for each subject
                                      cons = input_name,
                                      con_name = contrast,
                                      script = matlab_spmbatch_spec_est_level2),
        format = "file"
      ),
      values = crossing(contrast = contrast_names,
                        model_type = c("boxcar", "endspike")) %>% 
        mutate(target_name = syms(paste("level2", contrast, model_type, sep = ".")),
               input_name = syms(paste("con", contrast, model_type, sep = ".")),
               # it doesn't actually vary across the rows but tar_eval is not receiving vars set in the global env
               # so this is a quick way to get the task info into the call, I hope
               task_bids = paste0("task-", task))
    ),
    additional_targets
  )
  
}
