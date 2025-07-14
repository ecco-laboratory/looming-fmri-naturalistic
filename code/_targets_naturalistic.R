# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tidyverse)
library(rlang)

# Set target options:
tar_option_set(
  packages = c("withr",
               "matlabr",
               "tidymodels",
               "plsmod",
               "discrim",
               "tidyverse",
               "magrittr",
               "glue",
               "rlang",
               "qualtRics",
               "ggseg",
               "ggsegGlasser",
               "ragg",
               "cowplot"), # packages that your targets need to run
  error = "trim",
  controller = crew.cluster::crew_controller_slurm(
    workers = 8,
    seconds_idle = 30,
    options_cluster = crew.cluster::crew_options_slurm(
      verbose = TRUE,
      script_lines = "#SBATCH --account=default",
      log_output = "/home/%u/log/crew_log_%A.out",
      log_error = "/home/%u/log/crew_log_%A.err",
      memory_gigabytes_required = 8,
      cpus_per_task = 1,
      time_minutes = 1339,
      partition = "day-long"
    )
  )
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source(c("code/R/utils/",
             "code/R/make-stimlist.R",
             "code/R/call-activations.R",
             "code/R/define-targets-fmri.R",
             "code/R/parse-ann-preds.R",
             "code/R/parse-qualtrics.R",
             "code/R/parse-events.R",
             "code/R/parse-confounds.R",
             "code/R/parse-bold.R",
             "code/R/call-spm.R",
             "code/R/call-canlabtools.R",
             "code/R/utils/tidymodels-helpers.R",
             "code/R/proc-post-canlabtools.R",
             "code/R/proc-classify-controlled.R",
             "code/R/parse-parcels.R",
             "code/R/make-plots.R"))
# Regular source this script because it's not called by a target per se
# but is necessary for target construction
# source("code/R/set-study-defaults.R")

## other useful global vars ----

# format fMRI parameters from set-study-defaults.R
task_defaults_list <- jsonlite::read_json(here::here("task_defaults.json"), simplifyVector = TRUE) %>% 
  as_tibble() %>% 
  filter(task == "task-naturalistic") %>% 
  select(-task) %>% 
  as.list()

# target definition

## targets: scripts (sometimes in other project folders) ----

# the script targets are defined externally so this can be called in both targets pipelines
source("code/R/define-targets-scripts.R")

targets_scripts_naturalistic <- list(
  tar_target(
    name = matlab_spmbatch_contrast_level1,
    command = here::here("code", "matlab", "calc_contrasts_level1_naturalistic.m"),
    format = "file"
  )
)

## targets: looming stimuli of various kinds ----

folder_videos <- here::here("ignore", "stimuli", "videos", "naturalistic")
# just the video IDs for the rotting food videos
# since they aren't otherwise tagged (none of the food looms!)
# any food videos with other video IDs are tasty
gross_food <- sprintf("%04d_01.mp4", 49:58)

targets_stimuli <- list(
  # this will really only be used to generate placeholder activations
  # to patch into a full encoding model activation timecourse for the naturalistic task
  # I made it by running a fixation screenshot in iMovie for 1 second like a clown
  tar_target(name = video_fixation,
             command = here::here("ignore", "stimuli", "videos", "fixation.mp4"),
             format = "file"),
  tar_files(name = videos,
            command = list.files(folder_videos, pattern = ".mp4", full.names = TRUE)),
  # these don't really need to get shown to humans, we can show them the native fps ones
  # this is for the encoding models
  tar_target(name = videos_10fps,
             command = {
               videos
               resample_video_fps(in_path = folder_videos,
                                  out_path = here::here("ignore", "stimuli", "videos", "naturalistic_10fps"),
                                  script = py_resample_video_fps)
             },
             format = "file"),
  tar_target(name = metadata_videos,
             command = {
               out_path <- file.path(folder_videos, "metadata.csv")
               video_paths <- paste(videos, collapse = " ")
               command_args <- c(py_get_video_metadata,
                                 paste("-i", video_paths),
                                 paste("-o", out_path))
               
               run_python_target(command_args, out_path, conda_path)
             },
             format = "file"),
  tar_target(name = annotations_videos,
             command = file.path(folder_videos, "annotations.csv"),
             format = "file"),
  tar_target(name = ck2017_imagenet_categories,
             command = here::here("ignore", "data", "norm", "ck2017_imagenet_categories.csv"),
             format = "file"),
  tar_target(imagenet_category_labels,
             command = here::here("ignore", "imagenet_categories_synset.csv"),
             format = "file")
)

## targets: fMRI stimlist materials for naturalistic jumpy videos task ----

targets_stimlists <- list(
  tar_target(name = stims,
             command = read_csv(metadata_videos) %>% 
               select(video, height, width, duration) %>% 
               left_join(read_csv(annotations_videos) %>% 
                           mutate(video_id = paste0(video_id, ".mp4")) %>% 
                           rename(animal_type = animal), 
                         by = c("video" = "video_id")) %>% 
               filter(animal_type %in% c("dog", "frog", "spider", "cat", "food"), is.na(ignore)) %>% 
               select(-ignore)
  ),
  tar_map(
    values = tibble(order_num = sprintf("order%02d", 1:2)),
    tar_target(name = stimlist,
               command = make_stimlist_naturalistic(stims, order_num),
               format = "file")
  )
)

## targets: encoding model predictions (depend only on stimuli) ----

targets_encoding.models <- list(
  tar_target(name = activations.flynet_fixation,
             command = calc_flynet_activations(videos = video_fixation,
                                               out_path = here::here("ignore",
                                                                     "data",
                                                                     "encoding",
                                                                     "activations.flynet_fixation.csv"),
                                               script = py_calc_flynet_activations,
                                               weights = weights_flynet),
             format = "file"),
  tar_target(name = activations.flynet_raw,
             command = calc_flynet_activations(videos = videos_10fps,
                                               out_path = here::here("ignore",
                                                                     "data",
                                                                     "encoding",
                                                                     "activations.flynet_naturalistic.csv"),
                                               script = py_calc_flynet_activations,
                                               weights = weights_flynet),
             format = "file"),
  tar_target(name = activations.alexnet_raw,
             command = calc_alexnet_activations(videos = videos_10fps,
                                                out_path = here::here("ignore",
                                                                      "data",
                                                                      "encoding",
                                                                      "activations.alexnet_naturalistic.csv"),
                                                script = py_calc_alexnet_activations),
             format = "file"),
  tar_target(name = hitprobs.flynet_raw,
             command = calc_flynet_activations(videos = videos_10fps,
                                               out_path = here::here("ignore",
                                                                     "data",
                                                                     "encoding",
                                                                     "hitprobs.flynet_naturalistic.csv"),
                                               script = py_calc_flynet_activations,
                                               weights = weights_flynet,
                                               output_type = "hit_probs"),
             format = "file"),
  tar_target(name = activations.onoff_raw,
             command = {
               out_path <- here::here("ignore",
                                      "data",
                                      "encoding",
                                      "activations.onoff_naturalistic.csv")
               
               activations.flynet_raw %>% 
                 read_csv(name_repair = "unique") %>% 
                 select(frame, video) %>% 
                 # just one unit. just one activation. just do it.
                 mutate(`0` = 1) %>% 
                 write_csv(file = out_path)
               
               out_path
             },
             format = "file")
)

# targets: MRIQC type stuff across all subjects, not just manually included ----

# use dynamic branching here instead of static branching because we just want every available fmriprepped confounds file
# not just those for a pre-set group of subjects whose targets need to be named
# this won't tell us which subjects are in there, mind you!
# but it _should_ update itself if new files appear that match the list.files pattern
targets_qc <- list(
  tar_files(name = confounds_allsubs,
            command = get_all_raw_confounds(task = "naturalistic")),
  tar_target(name = signal_quality,
             command = get_labeled_noise_measures_by_subject(confounds_allsubs,
                                                             tr_duration = task_defaults_list$tr_duration,
                                                             disdaq_duration = task_defaults_list$disdaq_duration),
             pattern = map(confounds_allsubs)),
  tar_target(name = signal_quality_summary,
             command = signal_quality %>% 
               group_by(subj_num, run_num) %>% 
               # greater than half a voxel
               summarize(n_spikes = sum(framewise_displacement > 1.35)) %>% 
               summarize(n_spike_runs = sum(n_spikes > 0), n_spikes_total = sum(n_spikes)) %>% 
               # 2025-04-07: EXCLUDE SUBJECTS WHO HAVE A SPIKE IN ALL THREE RUNS
               # if you look, those subjects are also the ones who are moving a fuckton
               arrange(desc(n_spike_runs))),
  tar_target(name = signal_quality_by_stims,
             command = signal_quality %>% 
               group_by(subj_num, run_num) %>% 
               mutate(tr_num = 1:n()) %>% 
               filter(tr_num > 16) %>% 
               group_by(subj_num) %>% 
               mutate(tr_num = 1:n()) %>% 
               ungroup() %>% 
               # TODO: events.timecourse_all.subs only contains included subjects
               # stretch goal: make an alternate version of that target that instead crawls every available beh events csv, including from "excluded" subjects
               # just like how signal_quality gets it
               left_join(events.timecourse_all.subs, by = c("subj_num", "tr_num")) %>% 
               filter(video_id != "fixation") %>% 
               group_by(video_id, subj_num) %>% 
               summarize(max_fd = max(framewise_displacement), 
                         max_log.fd = max(log(framewise_displacement))) %>% 
               left_join(stims %>% 
                           select(video, has_loom, animal_type), 
                         by = c("video_id" = "video"))
             ),
  tar_target(name = plot_signal_quality_by_stims,
             command = ggplot(signal_quality_by_stims,
                              aes(x = animal_type, 
                                  y = max_fd, 
                                  color = factor(has_loom))) + 
               geom_hline(yintercept = 0, linetype = "dotted") + 
               geom_jitter(width = 0.1, alpha = 0.1) + 
               geom_pointrange(stat = "summary", fun.data = "mean_se"))
)

# targets: maps out by subject x run ----

# Use the participants.tsv file as a code-agnostic way of tracking which subjects to use
# It must be edited MANUALLY to label subjects as group "use" once their fmriqc has been checked and approved
# that way, only subjects manually approved will be included in these analyses
participants <- inject(here::here(!!!path_here_fmri, "participants.tsv")) %>% 
  read_tsv(comment = "#") %>% 
  # Anyone who didn't do all 3 runs should be marked as unusable for naturalistic
  filter(group %in% c("use_both", "use_naturalistic")) %>% 
  select(subject = participant_id)

# TODO 2025-04-03: can we decode object category directly from AlexNet activations? which category units activate in response to the different stimuli?
subtargets_encoding.timecourses_by.run <- list(
  # 2025-04-04: Phil wants to be able to look at these so output them to CSV instead
  tar_target(events.timecourse,
             command = make_condition_timecourse(onsets = events,
                                                 out_path = here::here("ignore", "outputs", sprintf("task-naturalistic_%s_%s_events.csv", subject, run_bids)),
                                                 tr_duration = task_defaults_list$tr_duration,
                                                 n_trs = task_defaults_list$n_trs_kept),
             format = "file"),
  tar_target(name = activations.flynet,
             command = make_encoding_timecourse_matlab(onsets = events,
                                                       path_stim_activations = activations.flynet_raw,
                                                       out_path = here::here("ignore", "outputs", sprintf("task-naturalistic_%s_%s_acts-flynet.csv", subject, run_bids)),
                                                       tr_duration = task_defaults_list$tr_duration,
                                                       run_duration = task_defaults_list$n_trs_kept * task_defaults_list$tr_duration),
             format = "file"),
  tar_target(name = activations.alexnet,
             command = make_encoding_timecourse_matlab(onsets = events,
                                                       path_stim_activations = activations.alexnet_raw,
                                                       out_path = here::here("ignore", "outputs", sprintf("task-naturalistic_%s_%s_acts-alexnet.csv", subject, run_bids)),
                                                       tr_duration = task_defaults_list$tr_duration,
                                                       run_duration = task_defaults_list$n_trs_kept * task_defaults_list$tr_duration),
             format = "file"),
  tar_target(name = activations.onoff,
             command = make_encoding_timecourse_matlab(onsets = events,
                                                       path_stim_activations = activations.onoff_raw,
                                                       out_path = here::here("ignore", "outputs", sprintf("task-naturalistic_%s_%s_acts-onoff.csv", subject, run_bids)),
                                                       tr_duration = task_defaults_list$tr_duration,
                                                       run_duration = task_defaults_list$n_trs_kept * task_defaults_list$tr_duration,
                                                       fixation_activation = tibble(`0` = 0)))
)

# attention! this is the innermost tar_map, which defines RUN-UNIQUE targets
targets_fmri_by.run <- make_targets_fmri_by.run(n_runs = task_defaults_list$n_runs,
                                                task = "naturalistic",
                                                additional_targets = subtargets_encoding.timecourses_by.run)

contrast_names <- c("dog",
                    "cat",
                    "frog",
                    "spider",
                    "food",
                    "looming",
                    "looming.baseline",
                    "stimuli")

# targets: maps out by task x subject (combines across run) ----

# we only need the full single-trial events dataframes for naturalistic
# to prepare to set the single-trial beta niftis as targets later
subtarget_events_combine <- tar_combine(name = events,
                                     targets_fmri_by.run[["events"]],
                                     command = bind_rows(!!!.x,
                                                         .id = "run") %>% 
                                       mutate(run = as.integer(run))
)

subtargets_encoding.timecourses_combine <- list(
  # 2025-04-04: the source is now output to CSV so this tar_combine needs to read_csv
  tar_combine(name = events.timecourse,
              targets_fmri_by.run[["events.timecourse"]],
              command = {
                vctrs::vec_c(!!!.x) %>% 
                map(\(x) read_csv(x)) %>% 
                bind_rows(.id = "target_name") %>% 
                  mutate(run_num = as.integer(str_sub(target_name, start = -2L)),
                         # bc this is eventually going to get bound with self-report data that has it just by video id
                         video_id = if_else(condition != "fixation",
                                            paste0(str_sub(condition, start = -7L), ".mp4"),
                                            condition),
                         tr_num = as.integer(tr_num + task_defaults_list$n_trs_kept * (run_num-1))) %>% 
                  select(-run_num, -target_name, -condition)
              }),
  tar_combine(name = activations.flynet,
              targets_fmri_by.run[["activations.flynet"]]),
  tar_combine(name = activations.alexnet,
              targets_fmri_by.run[["activations.alexnet"]]),
  tar_combine(name = activations.onoff,
              targets_fmri_by.run[["activations.onoff"]])
)

subtargets_bold.masked <- list(
  tar_target(name = bold.masked.sc,
             command = canlabtools_mask_fmri_data(out_path = here::here("ignore", "data", "canlabtools", 
                                                                        sprintf("%s_%s_region-sc_bold.mat", subject, task_bids)),
                                                  tr_duration = task_defaults_list$tr_duration,
                                                  trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                                  bolds = bold.smoothed,
                                                  confounds = confounds.prespm,
                                                  roi = "Bstem_SC",
                                                  script = matlab_mask_fmri_data),
             format = "file"),
  tar_target(name = bold.masked.amyg,
             command = canlabtools_mask_fmri_data(out_path = here::here("ignore", "data", "canlabtools", 
                                                                        sprintf("%s_%s_region-amyg_bold.mat", subject, task_bids)),
                                                  tr_duration = task_defaults_list$tr_duration,
                                                  trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                                  bolds = bold.smoothed,
                                                  confounds = confounds.prespm,
                                                  roi = "Amyg",
                                                  script = matlab_mask_fmri_data),
             format = "file")
)

targets_fmri_by.subject <- make_targets_fmri_by.subject(participants,
                                                        targets_fmri_by.run,
                                                        contrast_names,
                                                        task = "naturalistic",
                                                        additional_targets = c(subtarget_events_combine, 
                                                                               subtargets_encoding.timecourses_combine,
                                                                               subtargets_bold.masked))

## targets: maps out by task x contrast (combines across subject) ----

### SPM level 2 group analyses ----
# this eval factory is for aggregating stuff across all subs
subtargets_fmri_across.subject <- list(
  tar_combine(name = activations.flynet_all.subs,
              targets_fmri_by.subject[["activations.flynet"]],
              command = list(!!!.x)),
  tar_combine(name = activations.alexnet_all.subs,
              targets_fmri_by.subject[["activations.alexnet"]],
              command = list(!!!.x)),
  tar_combine(name = activations.onoff_all.subs,
              targets_fmri_by.subject[["activations.onoff"]],
              command = list(!!!.x)),
  tar_combine(name = bold.smoothed_all.subs,
              targets_fmri_by.subject[["bold.smoothed"]],
              command = list(!!!.x)),
  tar_combine(name = confounds.prespm_all.subs,
              targets_fmri_by.subject[["confounds.prespm"]],
              command = list(!!!.x)),
  # the components of these are already one per subject so we don't need to keep them as list
  tar_combine(name = bold.masked.sc_all.subs,
              targets_fmri_by.subject[["bold.masked.sc"]]),
  tar_combine(name = bold.masked.amyg_all.subs,
              targets_fmri_by.subject[["bold.masked.amyg"]]),
  tar_combine(name = events.timecourse_all.subs,
              targets_fmri_by.subject[["events.timecourse"]],
              command = bind_rows(!!!.x, .id = "target_name") %>% 
                mutate(subj_num = as.integer(str_sub(target_name, start = -4L))) %>% 
                select(-target_name) %>% 
                nest(events = -subj_num) %>% 
                mutate(fold_num = 1:n()) %>% 
                unnest(events)
              ),
  tar_target(events.timecourse.endspike_all.subs,
             command = events.timecourse_all.subs %>% 
               nest(events = -c(subj_num, fold_num)) %>% 
               mutate(events = map(events, \(x) relabel_timecourse_endspike(x))) %>% 
               unnest(events))
)

### canlabtools-based group analyses so help me god ----

subtargets_fmri_canlabtools_by.subject <- tar_map(
  values = participants %>% 
    mutate(bolds = sprintf("bold.smoothed_%s", subject),
           confounds = sprintf("confounds.prespm_%s", subject),
           across(c(bolds, confounds), \(x) syms(str_replace_all(x, "-", "\\.")))),
  # model-based connectivity is set up with single-subject targets to start
  # that depend on the all-subjects pred.encoding.xval target
  # to reduce memory load so that only one subject's whole-brain BOLD needs to be loaded in at a time to calculate the correlations
  tar_target(name = wb.model.connectivity,
             command = canlabtools_fit_model_connectivity(out_path = here::here("ignore", "outputs", sprintf("%s_naturalistic_wb.conn.%s_sc.csv", subject, encoding_type_full)),
                                                          tr_duration = task_defaults_list$tr_duration,
                                                          trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                                          bolds = bolds,
                                                          confounds = confounds,
                                                          pred.encoding.roi = pred.encoding.xval,
                                                          script = matlab_fit_model_connectivity),
             format = "file"),
  names = subject
)

subtargets_fmri_canlabtools_by.model <- tar_map(
  values = tribble(~encoding_type1, ~encoding_type2,
                   "flynet", "only",
                   "alexnet", "only",
                   "flynet", "alexnet",
                   "onoff", "only",
                   "flynet", "onoff",
                   "alexnet", "onoff") %>% 
    mutate(encoding_type_full = paste(encoding_type1, encoding_type2, sep = "."),
           activations1_all.subs = syms(sprintf("activations.%s_all.subs", encoding_type1)),
           activations2_all.subs = syms(sprintf("activations.%s_all.subs", encoding_type2)),
           activations2_all.subs = map_if(activations2_all.subs, 
                                          encoding_type2 == "only",
                                          \(x) NA)),
  tar_target(name = encoding.xval,
             command = canlabtools_fit_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.%s_sc.csv", encoding_type_full)),
                                                    out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.%s_sc.csv", encoding_type_full)),
                                                    out_path_betas = here::here("ignore", "outputs", sprintf("naturalistic_betas.%s_sc.csv", encoding_type_full)),
                                                    tr_duration = task_defaults_list$tr_duration,
                                                    bolds = bold.masked.sc_all.subs,
                                                    # by default fits one encoding model at a time
                                                    activations1 = activations1_all.subs,
                                                    activations2 = activations2_all.subs,
                                                    # important to do it as a list so they'll go into separate cells
                                                    # and also in the same alphabetical order as the out filenames!
                                                    script = matlab_fit_pls),
             format = "file"),
  tar_target(name = perf.encoding.xval,
             command = read_csv(encoding.xval[grepl("perf", encoding.xval)], col_names = FALSE)),
  tar_target(name = pred.encoding.xval,
             command = encoding.xval[grepl("pred", encoding.xval)],
             format = "file"),
  # TODO 2025-04-03: 1. Look at which AlexNet units are 
  tar_target(name = betas.encoding.xval,
             command = encoding.xval[grepl("betas", encoding.xval)],
             format = "file"),
  tar_target(name = betas.encoding.overall,
             command = canlabtools_fit_encoding_pls_no.xval(out_path_betas = here::here("ignore", "outputs", sprintf("naturalistic_betas.noxval.%s_sc.csv", encoding_type_full)),
                                                            tr_duration = task_defaults_list$tr_duration,
                                                            bolds = bold.masked.sc_all.subs,
                                                            # by default fits one encoding model at a time
                                                            activations1 = activations1_all.subs,
                                                            activations2 = activations2_all.subs,
                                                            # important to do it as a list so they'll go into separate cells
                                                            # and also in the same alphabetical order as the out filenames!
                                                            script = matlab_fit_pls_no_xval),
             format = "file"),
  tar_target(name = statmap.perf.encoding,
             command = {
               # canlabtools fmri_data write method appears to forcibly change periods in file names to underscores
               # so we have to feed in something that won't get changed so that the actual output path matches the expected one from here
               this_encoding_label <- str_replace(encoding_type_full, "\\.", "_")
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_perf_%s_sc.nii", this_encoding_label)),
                                          roi = "Bstem_SC",
                                          values = summarize_tvals_pre_statmap(perf.encoding.xval),
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_map(
    values = crossing(events_type = c("events.raw", "events.endspike"),
                      outcome_category = c("loom", "obj", "obj.loom")) %>% 
      mutate(events_label = case_match(events_type, "events.raw" ~ "", "events.endspike" ~ ".endspike"),
             events_target = syms(sprintf("events.timecourse%s_all.subs", events_label))),
    tar_target(name = encoding.decoding,
               command = pred.encoding.xval %>%
                 fit_object_by_pattern(events_allsubs = left_join(events_target,
                                                                  beh, 
                                                                  by = c("subj_num", "video_id")),
                                       n_trs_kept_per_run = task_defaults_list$n_trs_kept,
                                       pattern_type = "encoding",
                                       outcome_categories = outcome_category,
                                       n_pls_comp = 20)
    ),
    tar_rep(name = perm.acc_encoding.decoding,
            command = encoding.decoding %>% 
              select(subj_num, preds) %>% 
              unnest(preds) %>% 
              filter(animal_type != "food") %>% 
              permute_acc_object_by_pattern(n_perms = 200),
            batches = 50,
            reps = 1
    ),
    tar_rep(name = perm.acc_by.category_encoding.decoding,
            command = encoding.decoding %>% 
              select(subj_num, preds) %>% 
              unnest(preds) %>% 
              filter(animal_type != "food") %>% 
              permute_acc_object_by_pattern(n_perms = 200,
                                            acc_grouping_cols = animal_type),
            batches = 50,
            reps = 1
    ),
    tar_target(name = encoding.decoding.nosplit,
               command = pred.encoding.xval %>%
                 fit_object_by_pattern(events_allsubs = left_join(events_target,
                                                                  beh, 
                                                                  by = c("subj_num", "video_id")),
                                       n_trs_kept_per_run = task_defaults_list$n_trs_kept,
                                       pattern_type = "encoding",
                                       outcome_categories = outcome_category,
                                       n_pls_comp = 20,
                                       xval = FALSE)
    ),
    names = c(outcome_category, events_type)
  ),
  tar_target(name = encoding.selfreport_events.raw,
             command = fit_selfreport_by_pattern(pred.encoding.xval, 
                                                 events.timecourse_all.subs,
                                                 beh)
  ),
  tar_target(name = encoding.selfreport.nosplit,
             command = pred.encoding.xval %>% 
               get_ratings_by_encoding_space(left_join(events.timecourse.endspike_all.subs,
                                                       beh, 
                                                       by = c("subj_num", "video_id"))) %>% 
               filter(animal_type != "food") %>% 
               # SO!!! We need to make sure to avoid train-test leakage from the original encoding model to subsequent readout models.
               # we can do this by either:
               # fitting a single overall PLS regression for each of the outcome ratings (I honestly think this is fine. we don't need cross-validated generalization perf)
               # changing this PLS pipeline so that it gets fit along with the main encoding PLS train-test split (I would really prefer not to do this)
               pivot_longer(cols = starts_with("rating"),
                            names_to = "rating_type",
                            values_to = "rating",
                            names_prefix = "rating_") %>% 
               nest(.by = rating_type) %>% 
               mutate(model = map(data, \(x) {
                 n_voxels <- sum(grepl("voxel", names(x)))
                 fit_pls_single(in_data = x,
                                x_prefix = "voxel",
                                y_prefix = "rating",
                                # so that it will be equivalent to multiple regression
                                num_comp = n_voxels,
                                rm_x = TRUE) %>% 
                   pluck("model")
               })) %>% 
               mutate(coefs = map(model, \(x) extract_pls_workflow_coefs(x))) %>% 
               select(rating_type, coefs)
  ),
  # SEE ABOVE FOR THE DEFINITION OF THE INDIVIDUAL SUBJECT CONNECTIVITY TARGETS
  subtargets_fmri_canlabtools_by.subject,
  tar_combine(name = wb.model.connectivity,
              subtargets_fmri_canlabtools_by.subject[["wb.model.connectivity"]],
              # output as file in case Phil wants to read them in himself
              command = {
                out_path <- here::here("ignore", "outputs", sprintf("naturalistic_wb.conn.%s_sc.csv", encoding_type_full))
                
                connectivity <- vctrs::vec_c(!!!.x) %>% 
                  map(\(x) read_csv(x, col_names = FALSE)) %>% 
                  bind_cols(.name_repair = "unique")
                
                cat("transposing subjects back onto the rows\n")
                connectivity %<>%
                  t()
                
                cat("converting back to dataframe\n")
                connectivity %<>% 
                  as.data.frame()
                
                cat("writing tibble to file\n")
                vroom::vroom_write(connectivity, file = out_path, delim = ",", col_names = FALSE, progress = TRUE)
                
                out_path
                },
              format = "file"),
  # putting this in a sep target bc this takes like 8-10 min to calculate per statmap bc of all the voxels
  tar_target(name = tvals.wb.model.connectivity,
             command = wb.model.connectivity %>% 
               vroom::vroom(delim = ",", col_names = FALSE, col_types = c(.default = "d")) %>% 
               summarize_tvals_pre_statmap()
  ),
  tar_target(name = statmap.wb.model.connectivity,
             command = {
               # canlabtools fmri_data write method appears to forcibly change periods in file names to underscores
               # so we have to feed in something that won't get changed so that the actual output path matches the expected one from here
               this_encoding_label <- str_replace(encoding_type_full, "\\.", "_")
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_%s_sc.nii", this_encoding_label)),
                                          roi = NULL,
                                          values = tvals.wb.model.connectivity,
                                          threshold_p = .01,
                                          positive_only = TRUE,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_target(name = parcels.wb.connectivity,
             command = canlabtools_parcellate_avg(out_path = here::here("ignore", "outputs", sprintf("naturalistic_parcel_conn_%s_sc.csv", encoding_type_full)),
                                                  path_connectivity_allsubs_1 = wb.model.connectivity,
                                                  script = matlab_parcellate_avg),
             format = "file"),
  tar_map(
    values = tibble(sig_name = c("ceko2022", "kragel2015", "zhou2021"),
                    sig_folder = c("2021_Ceko_MPA2_multiaversive", "2015_Kragel_emotionClassificationBPLS", "2021_Zhou_Subjective_Fear")),
    tar_target(name = sig.sim,
               command = canlabtools_apply_wb_signature(out_path = here::here("ignore", "outputs", sprintf("naturalistic_sigsim-%s_%s_sc.csv", sig_name, encoding_type_full)),
                                                        fmri_data = wb.model.connectivity,
                                                        pattern_subdir = sig_folder,
                                                        script = matlab_apply_wb_signature),
               format = "file"),
    names = sig_name
  ),
  names = encoding_type_full
)

subtargets_fmri_canlabtools_compare.models <- list(
  subtargets_fmri_canlabtools_by.model,
  # this combined target is up here because it doesn't depend on any of the semi-recursive model comparison ones
  tar_combine(name = betas.encoding.xval,
              subtargets_fmri_canlabtools_by.model[["betas.encoding.xval"]],
              command = {
                vctrs::vec_c(!!!.x) %>% 
                  map(\(x) read_csv(x, col_names = FALSE) %>% 
                        # OMIT the intercept from the PLS
                        slice(-1)) %>% 
                  bind_rows(.id = "target_name") %>% 
                  rename(fold_num = X1) %>% 
                  mutate(encoding_type = str_split_i(target_name, "_", -1L)) %>% 
                  select(-target_name)
              }),
  #### decoding-from-encoding analyses mapped by decoding quantity and events labeling type ----
  tar_map(
    values = crossing(events_type = c("events.raw", "events.endspike"),
                      outcome_category = c("loom", "obj", "obj.loom")) %>% 
      mutate(events_label = case_match(events_type, "events.raw" ~ "", "events.endspike" ~ ".endspike"),
             events_target = syms(sprintf("events.timecourse%s_all.subs", events_label))),
    tar_target(name = bold.object,
               command = bold.masked.sc_all.subs %>%
                 fit_object_by_pattern(events_allsubs = left_join(events_target,
                                                                  beh, 
                                                                  by = c("subj_num", "video_id")),
                                       n_trs_kept_per_run = task_defaults_list$n_trs_kept,
                                       pattern_type = "bold",
                                       outcome_categories = outcome_category,
                                       n_pls_comp = 20)
    ),
    tar_rep(name = perm.acc_bold.object,
            command = bold.object %>% 
              select(preds) %>% 
              unnest(preds) %>% 
              permute_acc_object_by_pattern(n_perms = 100),
            batches = 10,
            reps = 1
    ),
    names = c(outcome_category, events_type)
  ),
  #### performance from joint encoding models with one set of activations knocked out ----
  tar_target(name = encoding.xval_flynet.minus.onoff,
             command = canlabtools_pred_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.flynet.minus.onoff_sc.csv")),
                                                     out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.flynet.minus.onoff_sc.csv")),
                                                     tr_duration = task_defaults_list$tr_duration,
                                                     bolds = bold.masked.sc_all.subs,
                                                     activations1 = activations.flynet_all.subs,
                                                     betas = betas.encoding.xval_flynet.onoff,
                                                     script = matlab_pred_pls),
             format = "file"),
  tar_target(name = encoding.xval_alexnet.minus.onoff,
             command = canlabtools_pred_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.alexnet.minus.onoff_sc.csv")),
                                                     out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.alexnet.minus.onoff_sc.csv")),
                                                     tr_duration = task_defaults_list$tr_duration,
                                                     bolds = bold.masked.sc_all.subs,
                                                     activations1 = activations.alexnet_all.subs,
                                                     betas = betas.encoding.xval_alexnet.onoff,
                                                     script = matlab_pred_pls),
             format = "file"),
  tar_target(name = perf.encoding.xval_flynet.minus.onoff,
             command = encoding.xval_flynet.minus.onoff %>% 
               .[grepl("perf", .)] %>% 
               read_csv(col_names = FALSE)),
  tar_target(name = pred.encoding.xval_flynet.minus.onoff,
             command = encoding.xval_flynet.minus.onoff %>% .[grepl("pred", .)],
             format = "file"),
  tar_target(name = perf.encoding.xval_alexnet.minus.onoff,
             command = encoding.xval_alexnet.minus.onoff %>% 
               .[grepl("perf", .)] %>% 
               read_csv(col_names = FALSE)),
  tar_target(name = pred.encoding.xval_alexnet.minus.onoff,
             command = encoding.xval_alexnet.minus.onoff %>% .[grepl("pred", .)],
             format = "file"),
  tar_target(name = statmap.perf.encoding_flynet.minus.alexnet,
             command = {
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_perf_flynet_minus_alexnet_sc.nii")),
                                          roi = roi_canlabtools,
                                          values = summarize_tvals_pre_statmap(perf.encoding.xval_flynet.only, perf.encoding.xval_alexnet.only),
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_target(name = parcels.wb.connectivity_flynet.minus.alexnet,
             command = canlabtools_parcellate_avg(out_path = here::here("ignore", "outputs", sprintf("naturalistic_parcel_conn_flynet_minus_alexnet_sc.csv")),
                                                  path_connectivity_allsubs_1 = wb.model.connectivity_flynet.only,
                                                  path_connectivity_allsubs_2 = wb.model.connectivity_alexnet.only,
                                                  fun_compare = "minus",
                                                  script = matlab_parcellate_avg),
             format = "file"),
  tar_target(name = tvals.wb.model.connectivity_flynet.minus.alexnet,
             command = {
               connectivity.flynet <- vroom::vroom(wb.model.connectivity_flynet.only, 
                                                   delim = ",", 
                                                   col_names = FALSE, 
                                                   col_types = c(.default = "d"))
               
               connectivity.alexnet <- vroom::vroom(wb.model.connectivity_alexnet.only, 
                                                    delim = ",", 
                                                    col_names = FALSE, 
                                                    col_types = c(.default = "d"))
               
               summarize_tvals_pre_statmap(in_data_1 = connectivity.flynet,
                                           in_data_2 = connectivity.alexnet,
                                           fun_compare = magrittr::subtract)
             }),
  tar_target(name = parcels.wb.connectivity_flynet.conj.alexnet,
             command = canlabtools_parcellate_avg(out_path = here::here("ignore", "outputs", sprintf("naturalistic_parcel_conn_flynet_conj_alexnet_sc.csv")),
                                                  path_connectivity_allsubs_1 = wb.model.connectivity_flynet.only,
                                                  path_connectivity_allsubs_2 = wb.model.connectivity_alexnet.only,
                                                  # this should work--matlab min(x, y) when x and y are arrays of identical dims
                                                  # appears to calculate the element-wise pmin which is what we want here
                                                  fun_compare = "min",
                                                  script = matlab_parcellate_avg),
             format = "file"),
  tar_target(name = tvals.wb.model.connectivity_flynet.conj.alexnet,
             command = {
               # 2025-04-09: Phil wants it this way (Tom Nichols citation?)
               # threshold the individual maps first
               tvals.flynet <- tvals.wb.model.connectivity_flynet.only %>% 
                 threshold_tvals_pre_statmap(threshold_p = .01)
               tvals.alexnet <- tvals.wb.model.connectivity_alexnet.only %>% 
                 threshold_tvals_pre_statmap(threshold_p = .01)
               
               pmin(tvals.flynet, tvals.alexnet)
               
             }),
  tar_target(name = statmap.wb.model.connectivity_flynet.minus.alexnet,
             command = {
               tvals_flynet <- vroom::vroom(wb.model.connectivity_flynet.only, 
                                            delim = ",", 
                                            col_names = FALSE, 
                                            col_types = c(.default = "d")) %>% 
                 summarize_tvals_pre_statmap() %>% 
                 threshold_tvals_pre_statmap(threshold_p = .01)
               
               tvals_alexnet <- vroom::vroom(wb.model.connectivity_alexnet.only, 
                                             delim = ",", 
                                             col_names = FALSE, 
                                             col_types = c(.default = "d")) %>% 
                 summarize_tvals_pre_statmap() %>% 
                 threshold_tvals_pre_statmap(threshold_p = .01)
               
               # there are many voxels where the difference is supra-threshold non-zero but the higher value isn't supra-threshold on its own
               # so for the diff statmap, zero out voxels ahead of time where the higher value on its own doesn't exceed 0
               # threshold_tvals_pre_statmap() retains only significant voxels as nonzero so this also has the effect 
               # of only keeping a - b diff voxels where a is significantly positive
               tvals_diff <- tvals.wb.model.connectivity_flynet.minus.alexnet
               tvals_diff[tvals_diff > 0 & tvals_flynet <= 0] <- 0
               tvals_diff[tvals_diff < 0 & tvals_alexnet <= 0] <- 0
               
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_flynet_minus_alexnet_sc.nii")),
                                          roi = NULL,
                                          values = tvals_diff,
                                          threshold_p = .01,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_target(name = statmap.wb.model.connectivity_flynet.conj.alexnet,
             command = {
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_flynet_conj_alexnet_sc.nii")),
                                          roi = NULL,
                                          values = tvals.wb.model.connectivity_flynet.conj.alexnet,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  # not a statmap! a preplot
  tar_target(name = wb.atlas_selected.rois,
             command = canlabtools_export_mask_nifti(out_path = here::here("ignore", "outputs", sprintf("mask_selected_rois.nii")),
                                                     # all of the Ctx_ ones to the Glasser atlas will match both left and right hemispheres
                                                     rois = list("Ctx_IFJ", # will match IFJa and IFJp
                                                                 "Ctx_FEF",
                                                                 c("Ctx_FFC", "Ctx_VVC"),
                                                                 c("Ctx_IPS1", "Ctx_MIP", "Ctx_LIP", "Ctx_VIP"),
                                                                 # the subcortical ones (other than SC which will come from Brainstem Navigator) are irreducibly bilateral
                                                                 "Amygdala",
                                                                 "Bstem_SC", 
                                                                 "Thal_Pulv",
                                                                 "Thal_LGN"),
                                                     script = matlab_export_mask_nifti),
             format = "file")
)

#### other downstream-of-canlabtools-encoding targets ----
subtargets_fmri_canlabtools_combined <- list(
  tar_combine(name = perf.encoding_combined,
              tar_select_targets(subtargets_fmri_canlabtools_compare.models, 
                                 starts_with("perf.encoding.xval")),
              command = {
                list(!!!.x) %>% 
                  # must average across voxels to consider data across ROIs
                  # but retain xval folds
                  map(\(x) x %>% 
                        rowwise() %>% 
                        mutate(perf = mean(c_across(everything()))) %>% 
                        ungroup() %>% 
                        mutate(fold_num = 1:n()) %>% 
                        select(fold_num, perf)) %>% 
                  bind_rows(.id = "target_name") %>% 
                  separate_wider_delim(cols = target_name, delim = "_", names = c(NA, "encoding_type"))
              }),
  tar_combine(name = pred.encoding.xval,
              tar_select_targets(subtargets_fmri_canlabtools_compare.models, 
                                 starts_with("pred.encoding.xval")),
              command = {
                vctrs::vec_c(!!!.x) %>% 
                  map(\(x) load_encoding_pred_allsubs(x)) %>% 
                  bind_rows(.id = "target_name") %>% 
                  mutate(encoding_type = str_split_i(target_name, "_", -2L)) %>% 
                  select(-target_name)
              }),
  tar_target(name = pcor.encoding_flynet.alexnet,
             command = calc_pcor_encoding_on_bold(bold.masked.sc_all.subs,
                                                  pred.encoding.xval_flynet.only,
                                                  pred.encoding.xval_alexnet.only,
                                                  "flynet",
                                                  "alexnet")
  ),
  tar_target(name = predata_auroc.cross_encoding.discrim.looming,
             command = encoding.decoding_loom_events.raw_flynet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               preproc_auroc_obj_from_loom()
  ),
  tar_target(name = auroc.cross_encoding.discrim.looming,
             command = predata_auroc.cross_encoding.discrim.looming %>% 
               roc_auc(truth = animal_value, .pred_outcome_loom)
  ),
  tar_rep(name = perm_auroc.cross_encoding.discrim.looming,
          command = encoding.decoding_loom_events.raw_flynet.only %>% 
            select(-coefs) %>% 
            unnest(preds) %>% 
            filter(animal_type != "food") %>% 
            permutations(permute = animal_type, times = 500) %>% 
            mutate(perm_auroc = map(splits, \(x) analysis(x) %>% 
                                      preproc_auroc_obj_from_loom() %>% 
                                      roc_auc(truth = animal_value, .pred_outcome_loom),
                                    .progress = "permuting")) %>% 
            select(-splits) %>% 
            unnest(perm_auroc),
          batches = 20,
          reps = 1
  ),
  tar_target(name = auroc.cross.allway_encoding.discrim.object,
             command = encoding.decoding_obj_events.raw_alexnet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               preproc_auroc_loom_from_obj() %>% 
               roc_auc(truth = has_loom, .pred_animal)
  ),
  tar_target(name = predata_auroc.cross_encoding.discrim.object,
             command = encoding.decoding_obj_events.raw_alexnet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               preproc_auroc_loom_from_obj() %>% 
               filter(animal_to_use != "food") %>% 
               group_by(animal_to_use)
  ),
  tar_target(name = auroc.cross_encoding.discrim.object,
             command = predata_auroc.cross_encoding.discrim.object %>% 
               roc_auc(truth = has_loom, .pred_animal)
  ),
  tar_rep(name = perm_auroc.cross_encoding.discrim.object,
          command = encoding.decoding_obj_events.raw_alexnet.only %>% 
            select(-coefs) %>% 
            unnest(preds) %>% 
            filter(animal_type != "food") %>% 
            permutations(permute = has_loom, times = 500) %>% 
            mutate(perm_auroc = map(splits, \(x) analysis(x) %>% 
                                      preproc_auroc_loom_from_obj() %>% 
                                      group_by(animal_to_use) %>% 
                                      roc_auc(truth = has_loom, .pred_animal),
                                    .progress = "permuting")) %>% 
            select(-splits) %>% 
            unnest(perm_auroc),
          batches = 20,
          reps = 1
  )
)

## targets imported from the controlled task ----
# these do NOT auto update when new subjects are processed for naturalistic
# you must tar_make the original controlled targets for that tracking
#  but if you do tar_make controlled, then the correct naturalistic dependency targets will be marked as out of date

store_controlled <- here::here("ignore", "_targets", "controlled")
metadata_videos_nback <- tar_read(metadata_videos_nback, store = store_controlled)
plot_controlled_norms <- tar_read(plot_norms, store = store_controlled)

targets_controlled <- list(
  tar_target(activations.flynet_raw_controlled,
             command = here::here("ignore",
                                  "data",
                                  "encoding",
                                  "activations.flynet_controlled.csv"),
             format = "file"),
  tar_target(activations.alexnet_raw_controlled,
             command = here::here("ignore",
                                  "data",
                                  "encoding",
                                  "activations.alexnet_controlled.csv"),
             format = "file"),
  tar_files(con.files_level1.5_controlled,
            command = list.files(here::here("ignore", "data", "canlabtools"), 
                                 pattern = "task-controlled_region-sc_con-.*ing\\..*\\.baseline.csv",
                                 full.names = TRUE)),
  tar_target(cons_level1.5_controlled,
             command = {
               file_path <- con.files_level1.5_controlled
               condition <- file_path %>% 
                 basename() %>% 
                 str_split_i("-", -1) %>% 
                 str_remove(".baseline.csv")
               
               read_csv(file_path, col_names = FALSE) %>% 
                 mutate(condition = condition) %>% 
                 separate_wider_delim(cols = condition, delim = ".", names = c("direction", "animal_type")) %>% 
                 select(direction, animal_type, subj_num = X1, everything())
               },
             pattern = map(con.files_level1.5_controlled)),
  tar_target(bold.object.controlled,
             command = classify_controlled(cons_level1.5_controlled, outcome_var = animal_type)),
  tar_target(bold.looming.controlled,
             command = classify_controlled(cons_level1.5_controlled, outcome_var = direction)),
  tar_target(perm_bold.object.controlled,
             command = permute_xval_classification(bold.object.controlled, n_perms = 10000)),
  tar_target(perm_bold.object.controlled_by.category,
             command = permute_xval_classification(bold.object.controlled, n_perms = 10000, grouping_cols = .obs)),
  tar_target(perm_bold.looming.controlled,
             command = permute_xval_classification(bold.looming.controlled, n_perms = 10000)),
  tar_target(summary_bold.object.controlled,
             command = calc_xval_perf(bold.object.controlled, df_perms = perm_bold.object.controlled)),
  tar_target(summary_bold.object.controlled_by.category,
             command = calc_xval_perf(bold.object.controlled, grouping_cols = .obs, df_perms = perm_bold.object.controlled_by.category)),
  tar_target(summary_bold.looming.controlled,
             command = calc_xval_perf(bold.looming.controlled, df_perms = perm_bold.looming.controlled)),
  tar_target(pattern.activation.controlled_flynet,
             command = calc_controlled_pattern_expression(cons_level1.5_controlled,
                                                          activations.flynet_raw_controlled,
                                                          metadata_videos_nback,
                                                          betas.encoding.xval_flynet.only)),
  tar_target(pattern.discrim.controlled_alexnet,
             command = calc_controlled_pattern_discrim(cons_level1.5_controlled,
                                                       activations.alexnet_raw_controlled,
                                                       metadata_videos_nback,
                                                       betas.encoding.xval_alexnet.only,
                                                       encoding.decoding_obj_events.endspike_alexnet.only))
)

## overall summary stat targets ----

summary_funs_bootstrap <- list(mean = \(x) mean(x, na.rm = TRUE),
                               sd = \(x) sd(x, na.rm = TRUE),
                               ci95.lower = \(x) quantile(x, .025, na.rm = TRUE),
                               ci95.upper = \(x) quantile(x, .975, na.rm = TRUE))

subtargets_fmri_summary <- list(
  tar_map(
    values = tibble(rating_type = c("pleasantness", "arousal", "fear")),
    tar_target(name = anova_beh,
               command = beh %>% 
                 filter(animal_type != "food") %>% 
                 mutate(animal_type = factor(animal_type, levels = c("dog", "cat", "frog", "spider")),
                        # effectively order and grand mean center it. 
                        # this also reduces the number of parameters to fit by not dummy-coding them separately
                        animal_type_num = as.numeric(animal_type) - 2.5) %>% 
                 rename(this_rating = paste0("rating_", rating_type)) %>% 
                 # we can now keep it MAXIMAL and the fit isn't singular bc fewer parameters to be collinear lol
                 lmerTest::lmer(this_rating ~ has_loom * animal_type_num + (1 + has_loom * animal_type_num | subj_num), data = .) %>% 
                 anova(ddf = "Kenward-Roger"))
  ),
  tar_target(name = summary_perf.encoding,
             command = perf.encoding_combined %>% 
               group_by(encoding_type) %>% 
               mutate(r2 = perf^2) %>% 
               rename(r = perf) %>% 
               summarize(across(c(r, r2), summary_funs_bootstrap)) %>% 
               mutate(r_cohens.d = cohens.d.2(r_mean, r_sd, 0, r_sd),
                      r2_cohens.d = cohens.d.2(r2_mean, r2_sd, 0, r2_sd))
  ),
  tar_target(name = summary_perf.encoding_alexnet.minus.flynet,
             command = perf.encoding_combined %>% 
               filter(encoding_type %in% c("flynet.only", "alexnet.only")) %>% 
               pivot_wider(names_from = encoding_type, values_from = perf) %>% 
               # a priori we know alexnet is higher so this will get it signed positive. whatever
               mutate(perf_diff = alexnet.only - flynet.only,
                      r2_diff = alexnet.only^2 - flynet.only^2) %>% 
               summarize(diff_mean = mean(perf_diff),
                         diff_se = sd(perf_diff)/sqrt(length(perf_diff)),
                         diff_cohens.d = cohens.d(alexnet.only, flynet.only),
                         r2.diff_mean = mean(r2_diff),
                         r2.diff_cohens.d = cohens.d(alexnet.only^2, flynet.only^2))
  ),
  tar_target(name = summary_pcor.encoding_flynet.alexnet,
             command = pcor.encoding_flynet.alexnet %>% 
               rename(pcor = estimate, encoding_type = term) %>% 
               group_by(encoding_type) %>% 
               summarize(across(pcor, c(summary_funs_bootstrap))) %>% 
               mutate(pcor_cohens.d = cohens.d.2(pcor_mean, pcor_sd, 0, pcor_sd))
             ),
  tar_target(name = encoding.selfreport_alexnet.flynet.cbound,
             command = bind_rows(alexnet = encoding.selfreport_events.raw_alexnet.only,
                                 flynet = encoding.selfreport_events.raw_flynet.only,
                                 .id = "model_type") %>% 
               select(subj_num, rating_type, model_type, preds) %>% 
               unnest(preds) %>% 
               # Clearly some code generalization failure thing is making this column of all NAs pop up. idk man
               select(-rating)),
  tar_rep(name = boots_encoding.selfreport,
          command = encoding.selfreport_alexnet.flynet.cbound %>% 
            nest(.by = c(model_type, rating_type)) %>% 
            mutate(data = map(data, \(x) group_bootstraps(x, group = subj_num, times = 2000))) %>%
            unnest(data) %>% 
            mutate(cor.overall = map_dbl(splits, \(x) x %>% 
                                           analysis() %$% 
                                           cor(.obs, .pred, method = "spearman"),
                                         .progress = "bootstrapping pred-obs correlations"),
                   cor.animal = map_dbl(splits, \(x) x %>% 
                                          analysis() %>% 
                                          filter(animal_type != "food") %$%
                                          cor(.obs, .pred, method = "spearman"),
                                        .progress = "bootstrapping pred-obs correlations"),
                   cor.food = map_dbl(splits, \(x) x %>% 
                                        analysis() %>% 
                                        filter(animal_type == "food") %$%
                                        cor(.obs, .pred, method = "spearman"),
                                      .progress = "bootstrapping pred-obs correlations FOOD ONLY")) %>% 
            select(-splits),
          batches = 5,
          reps = 1
  ),
  tar_target(name = summary_encoding.selfreport,
             command = boots_encoding.selfreport %>%
               pivot_longer(cols = starts_with("cor."), 
                            names_to = "data_subset", 
                            values_to = "correlation", 
                            names_prefix = "cor.") %>% 
               group_by(data_subset, rating_type, model_type) %>% 
               summarize(across(correlation, 
                                summary_funs_bootstrap),
                         .groups = "drop")
  ),
  tar_rep(name = boots_encoding.selfreport_pcor,
          command = encoding.selfreport_alexnet.flynet.cbound %>% 
            pivot_wider(names_from = model_type, values_from = .pred) %>% 
            nest(.by = rating_type) %>% 
            mutate(data = map(data, \(x) group_bootstraps(x, group = subj_num, times = 2000))) %>%
            unnest(data) %>% 
            mutate(coefs.overall = map(splits, \(x) x %>% 
                                         analysis() %>% 
                                         mutate(alexnet = dense_rank(alexnet),
                                                flynet = dense_rank(flynet)) %>% 
                                         lm(scale(.obs) ~ scale(alexnet) + scale(flynet), data = .) %>% 
                                         broom::tidy(),
                                       .progress = "bootstrapping partial correlations"),
                   coefs.animal = map(splits, \(x) x %>% 
                                        analysis() %>% 
                                        # 2025-04-05: SPEARMAN CORRELATION because the split-half models are (predictably) way noisier out of sample. so there are some wild outliers
                                        mutate(alexnet = dense_rank(alexnet),
                                               flynet = dense_rank(flynet)) %>% 
                                        filter(animal_type != "food") %>% 
                                        # to estimate les partial correlations
                                        lm(scale(.obs) ~ scale(alexnet) + scale(flynet), data = .) %>% 
                                        broom::tidy(),
                                      .progress = "bootstrapping partial correlations"),
                   coefs.food = map(splits, \(x) x %>% 
                                      analysis() %>% 
                                      mutate(alexnet = dense_rank(alexnet),
                                             flynet = dense_rank(flynet)) %>% 
                                      filter(animal_type == "food") %>% 
                                      lm(scale(.obs) ~ scale(alexnet) + scale(flynet), data = .) %>% 
                                      broom::tidy(),
                                    .progress = "bootstrapping partial correlations FOOD ONLY")) %>% 
            select(-splits),
          batches = 5,
          reps = 1
  ),
  tar_target(name = summary_encoding.selfreport_pcor,
             command = boots_encoding.selfreport_pcor %>% 
               pivot_longer(cols = starts_with("coefs"), names_to = "data_subset", values_to = "coefs", names_prefix = "coefs.") %>% 
               unnest(coefs) %>% 
               filter(term != "(Intercept)") %>% 
               mutate(model_type = str_sub(term, start = 7L, end = -2L)) %>% 
               group_by(rating_type, data_subset, model_type) %>% 
               rename(pcor = estimate) %>% 
               summarize(across(pcor, summary_funs_bootstrap),
                         .groups = "drop")
               ),
  tar_target(name = summary_encoding.discrim.looming,
             command = calc_perm_pval_object_by_pattern(preds = encoding.decoding_loom_events.raw_flynet.only, 
                                                        perms = perm.acc_encoding.decoding_loom_events.raw_flynet.only)
  ),
  tar_target(name = summary_encoding.discrim.looming_by.category,
             command = calc_perm_pval_object_by_pattern(preds = encoding.decoding_loom_events.raw_flynet.only, 
                                                        perms = perm.acc_by.category_encoding.decoding_loom_events.raw_flynet.only,
                                                        grouping_cols = animal_type)
  ),
  tar_target(name = summary_encoding.discrim.looming_diff.spider,
             command = {
               acc_true <- encoding.decoding_loom_events.raw_flynet.only %>% 
                 select(-coefs) %>% 
                 unnest(preds) %>% 
                 calc_metrics_nested_classprobs(grouping_cols = animal_type) %>% 
                 filter(.metric == "roc_auc")
               
               auroc_diff_spider <- acc_true$.estimate[acc_true$animal_type == "spider"] - mean(acc_true$.estimate[acc_true$animal_type != "spider"])
               
               acc_perm <- perm.acc_by.category_encoding.decoding_loom_events.raw_flynet.only %>% 
                 mutate(auroc_diff_spider_perm = map_dbl(acc_perm, \(x) {
                   x_filtered <- filter(x, .metric == "roc_auc")
                   x_filtered$.estimate[x_filtered$animal_type == "spider"] - mean(x_filtered$.estimate[x_filtered$animal_type != "spider"])
                 }),
                 auroc_diff_spider_true = auroc_diff_spider) %>% 
                 summarize(estimate_real = unique(auroc_diff_spider_true),
                           # need this to report permuted chance levels
                           estimate_perm = median(auroc_diff_spider_perm),
                           pval = (sum(auroc_diff_spider_perm > auroc_diff_spider_true)+1)/(n()+1))
               
               acc_perm
             }
  ),
  tar_target(name = summary_auroc.cross_encoding.discrim.looming,
             command = left_join(perm_auroc.cross_encoding.discrim.looming, auroc.cross_encoding.discrim.looming, 
                                 by = c("animal_to_guess", ".metric", ".estimator"),
                                 suffix = c("_perm", "_real")) %>% 
               group_by(animal_to_guess) %>% 
               summarize(pval = (sum(.estimate_perm > .estimate_real)+1)/(n()+1),
                         .estimate_real = unique(.estimate_real),
                         # need this to report permuted chance levels
                         .estimate_perm = median(.estimate_perm)) %>% 
               select(animal_to_guess, .estimate_real, .estimate_perm, pval)),
  tar_target(name = summary_encoding.discrim.object,
             command = calc_perm_pval_object_by_pattern(preds = encoding.decoding_obj_events.raw_alexnet.only %>% 
                                                          mutate(preds = map(preds, \(x) filter(x, animal_type != "food"))), 
                                                        perms = perm.acc_encoding.decoding_obj_events.raw_alexnet.only)
  ),
  tar_target(name = summary_auroc.cross_encoding.discrim.object,
             command = left_join(perm_auroc.cross_encoding.discrim.object, auroc.cross_encoding.discrim.object, 
                                 by = c("animal_to_use", ".metric", ".estimator"),
                                 suffix = c("_perm", "_real")) %>% 
               group_by(animal_to_use) %>% 
               summarize(pval = (sum(.estimate_perm > .estimate_real)+1)/(n()+1),
                         .estimate_real = unique(.estimate_real),
                         # need this to report permuted chance levels
                         .estimate_perm = median(.estimate_perm)) %>% 
               select(animal_to_use, .estimate_real, .estimate_perm, pval)),
  tar_target(name = summary_cons_level1.5_controlled,
             command = cons_level1.5_controlled %>% 
               rowwise() %>% 
               # average across voxels within subject x beta map
               mutate(con_mean = mean(c_across(starts_with("X")))) %>% 
               ungroup() %>% 
               select(-starts_with("X")) %>% 
               group_by(direction, subj_num) %>% 
               # average across animal types within subject
               summarize(con_mean_2 = mean(con_mean)) %>%
               pivot_wider(names_from = direction, values_from = con_mean_2) %>% 
               mutate(diff = looming - receding) %>% 
               summarize(across(c(looming, receding, diff), .fns = summary_funs_bootstrap),
                         diff_cohens.d = cohens.d(looming, receding))
             ),
  tar_target(name = summary_pattern.activation.controlled_flynet,
             command = pattern.activation.controlled_flynet %>% 
               select(where(\(x) !is.list(x))) %>%
               group_by(direction, subj_num) %>% 
               # average across animal types within subject
               summarize(mean_expression = mean(pattern_expression)) %>%
               pivot_wider(names_from = direction, values_from = mean_expression) %>% 
               mutate(diff = looming - receding) %>% 
               summarize(across(c(looming, receding, diff), .fns = summary_funs_bootstrap),
                         diff_cohens.d = cohens.d(looming, receding))),
  tar_target(name = summary_pattern.activation.controlled_by.object_flynet,
             command = pattern.activation.controlled_flynet %>% 
               select(where(\(x) !is.list(x))) %>%
               mutate(is_spider = if_else(animal_type == "spider", "spider", "other")) %>% 
               group_by(is_spider, direction, subj_num) %>% 
               summarize(mean_expression = mean(pattern_expression)) %>%
               pivot_wider(names_from = direction, values_from = mean_expression) %>% 
               mutate(diff = looming - receding) %>% 
               select(-looming, -receding) %>% 
               pivot_wider(names_from = is_spider, values_from = diff) %>% 
               mutate(diff = spider - other) %>% 
               summarize(across(c(spider, other, diff), .fns = summary_funs_bootstrap),
                         diff_cohens.d = cohens.d(spider, other))),
  tar_target(generalization.controlled_flynet,
             command = calc_controlled_perf(pattern.activation.controlled_flynet)),
  tar_target(name = wb.connectivity_sig.sim,
             command = {
               ceko2022 <- bind_rows(looming = read_csv(sig.sim_ceko2022_flynet.only),
                                     object = read_csv(sig.sim_ceko2022_alexnet.only),
                                     .id = "model_type") %>% 
                 rename_with(.fn = \(x) paste0("aversive_", tolower(str_split_i(x, "_", 1))), .cols = -model_type) %>% 
                 group_by(model_type) %>% 
                 mutate(fold_num = 1:n()) %>% 
                 ungroup()
               
               kragel2015 <- bind_rows(looming = read_csv(sig.sim_kragel2015_flynet.only),
                                       object = read_csv(sig.sim_kragel2015_alexnet.only),
                                       .id = "model_type") %>% 
                 rename_with(.fn = \(x) str_split_i(x, "_", 3), .cols = -model_type) %>% 
                 rename(amusement = amused, anger = angry, contentment = content, fear = fearful, sadness = sad, surprise = surprised) %>% 
                 rename_with(.fn = \(x) paste0("emotion_", x), .cols = -model_type) %>% 
                 group_by(model_type) %>% 
                 mutate(fold_num = 1:n()) %>% 
                 ungroup()
               
               zhou2021 <- bind_rows(looming = read_csv(sig.sim_zhou2021_flynet.only),
                                     object = read_csv(sig.sim_zhou2021_alexnet.only),
                                     .id = "model_type") %>% 
                 # it's a subjective fear signature but don't want the name to conflict with fear from kragel 2015
                 rename(VIFS_fear = VIFS.nii) %>% 
                 group_by(model_type) %>% 
                 mutate(fold_num = 1:n()) %>% 
                 ungroup()
               
               reduce(list(ceko2022, kragel2015, zhou2021), 
                      \(x, y) full_join(x, y, by = c("model_type", "fold_num")))
             }),
  tar_target(name = summary_wb.connectivity_sig.sim,
             command = wb.connectivity_sig.sim %>% 
               group_by(model_type) %>% 
               # there are still NaNs even though I told canlabtools to ignore missing before this...
               summarize(across(-fold_num, list(mean = \(x) mean(x, na.rm = TRUE), 
                                                sd = \(x) sd(x, na.rm = TRUE), 
                                                se = \(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.nan(x))))))
  ),
  tar_target(name = beta.comparison_category.selfreport,
             command = compare_betas_encoding_category_selfreport(encoding.decoding.nosplit_loom_events.raw_flynet.only,
                                                                  encoding.decoding.nosplit_obj_events.raw_alexnet.only,
                                                                  encoding.selfreport.nosplit_flynet.only,
                                                                  encoding.selfreport.nosplit_alexnet.only,
                                                                  n_bootstraps = 10000))
)

targets_fmri_across.subject <- make_targets_fmri_across.subject(targets_fmri_by.subject,
                                                                contrast_names,
                                                                task = "naturalistic",
                                                                additional_targets = c(subtargets_fmri_across.subject,
                                                                                       subtargets_fmri_canlabtools_compare.models,
                                                                                       subtargets_fmri_canlabtools_combined,
                                                                                       subtargets_fmri_summary))

## targets: behavioral data ----

targets_beh <- list(
  tar_target(name = norms_qualtrics_all.stimuli,
             command = get_splat_stimulus_norms_qualtrics()),
  tar_target(name = norms_qualtrics,
             command = join_raw_norms_to_stim_labels(norms_qualtrics_all.stimuli, 
                                                     read_csv(annotations_videos) %>% 
                                                       select(-ignore) %>% 
                                                       rename(animal_type = animal),
                                                     loom_col = "has_loom")
  ),
  tar_files(name = beh.raw,
            # now filters for files belonging to subjects marked as usable in participants.tsv
            command = {
              all_files <- list.files(here::here("ignore", "data", "beh"),
                                 pattern = "task-naturalistic_beh", 
                                 full.names = TRUE, 
                                 recursive = TRUE)
              
              all_subjects <- str_split_i(all_files, "/", 9)
              all_files[all_subjects %in% participants$subject]
              }),
  # these two targets map over every beh.raw file so you will need to filter for kept subjects later
  tar_target(name = beh,
             command = beh.raw %>% 
               read_csv() %>% 
               select(subj_num = participant, video_id, has_loom, animal_type, ends_with("rating")) %>% 
               # turn the suffix into a prefix
               rename_with(\(x) str_split_i(x, "_", 1) %>% 
                             paste0("rating_", .),
                           ends_with("rating")) %>% 
               mutate(subj_num = as.integer(subj_num)) %>% 
               filter(!is.na(video_id)),
             pattern = map(beh.raw)),
  # 2025-02-21: currently returning whole tidy RDM (not halved in case you want to reorder when plotting, like with hclust)
  tar_target(name = rdm.beh,
             command = beh %>% 
               select(-has_loom, -animal_type) %>% 
               expand(subj_num,
                      nesting(video1 = video_id, 
                              pleasantness1 = rating_pleasantness, 
                              arousal1 = rating_arousal, 
                              fear1 = rating_fear), 
                      nesting(video2 = video_id, 
                              pleasantness2 = rating_pleasantness, 
                              arousal2 = rating_arousal, 
                              fear2 = rating_fear)) %>% 
               filter(video1 != video2) %>% 
               # halve_tidy_rdm(row_col = video1, col_col = video2) %>% 
               mutate(diff_pleasantness = abs(pleasantness1 - pleasantness2),
                      diff_arousal = abs(arousal1 - arousal2),
                      diff_fear = abs(fear1 - fear2)) %>% 
               select(subj_num, video1, video2, starts_with("diff")),
             pattern = map(beh)
  )
)

## targets: plots for showing ----

# Set a font that should come pre-installed with Linux on all the compute nodes
# Been sort of struggling to get Google Fonts to show up even with the agg_png write device
# But this ought to look close enough to whatever other Microsoft fonts are used in figures
this_theme <- theme_bw(base_size = 16, base_family = "Nimbus Sans") +
  theme(plot.background = element_blank(),
        legend.background = element_blank())

theme_brain <- theme_void(base_size = 16) +
  theme(plot.background = element_blank(),
        legend.background = element_blank())

targets_plots <- list(
  tar_target(plot_ratings,
             command = beh %>% 
               relabel_cols_for_plot_naturalistic() %>% 
               plot_selfreport_ratings() +
                 this_theme
  ),
  tar_target(plot_ratings_nofood,
             command = beh %>% 
               filter(animal_type != "food") %>% 
               relabel_cols_for_plot_naturalistic() %>% 
               plot_selfreport_ratings() +
               this_theme
  ),
  tar_target(plot_schematic_unit.activation,
             # since it's a schematic, just use subject 1
             command = list("looming" = activations.flynet_sub.0001,
                            "object" = activations.alexnet_sub.0001) %>% 
               plot_sample_timecourse_activation() +
               this_theme +
               theme(aspect.ratio = 1/3)
             ),
  tar_target(plot_schematic_encoding.bold.pred,
             command = list("looming" = pred.encoding.xval_flynet.only,
                            "object" = pred.encoding.xval_alexnet.only) %>% 
               plot_sample_timecourse_encoding() +
               this_theme +
               # the aspect ratio is for a single panel on faceted plots
               theme(aspect.ratio = 1/3)
  ),
  tar_target(plot_schematic_bold,
             command = plot_sample_timecourse_bold(bold.masked.sc_sub.0001) +
               this_theme +
               theme(aspect.ratio = 1/3)
  ),
  tar_target(plot_bold.object_auroc,
             command = plot_auroc_8cat(bold.object_8cat_events.endspike) +
               this_theme
  ),
  tar_target(plot_encoding.perf,
             command = plot_encoding_performance(perf.encoding_combined) +
               this_theme
  ),
  tar_target(plot_encoding.perf_flynet.alexnet,
             command = plot_encoding_performance(perf.encoding_combined, 
                                                 encoding_types = c("looming" = "flynet.only", 
                                                                    "object" = "alexnet.only")) +
               guides(color = "none") + 
               labs(x = NULL,
                    subtitle = NULL) +
               this_theme +
               theme(axis.title = element_text(size = rel(1.3)),
                     axis.text = element_text(size = rel(1.05)))
  ),
  tar_target(plot_encoding.pcor_flynet.alexnet,
             command = plot_encoding_performance_pcor(pcor.encoding_flynet.alexnet) +
               this_theme
  ),
  tar_target(plot_encoding.confusion.looming_flynet,
             command = encoding.decoding_loom_events.raw_flynet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               mutate(across(c(.obs, .pred), \(x) fct_recode(x,
                                                             "looming" = "loom",
                                                             "no looming" = "no.loom"))) %>% 
               plot_confusion() +
               this_theme +
               theme(aspect.ratio = 1)
  ),
  tar_target(plot_encoding.confusion.looming_flynet_events.raw,
             command = encoding.decoding_loom_events.raw_flynet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               mutate(across(c(.obs, .pred), \(x) fct_recode(x,
                                                             "looming" = "loom",
                                                             "no looming" = "no.loom"))) %>% 
               plot_confusion() +
               this_theme +
               theme(aspect.ratio = 1)
  ),
  tar_target(plot_encoding.confusion.object,
             command = bind_rows(looming = encoding.decoding_obj_events.raw_flynet.only,
                                 object = encoding.decoding_obj_events.raw_alexnet.only,
                                 .id = "model_type") %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               mutate(across(c(.obs, .pred), \(x) x %>% 
                               fct_relevel("dog", "cat", "frog", "spider")),
                      has_loom = if_else(has_loom == 1, "looming", "no looming"),
                      model_type = paste0(model_type, " encoding")) %>% 
               plot_confusion(facet_var_row = has_loom,
                                   facet_var_col = model_type) +
               this_theme +
               theme(aspect.ratio = 1)
  ),
  tar_target(plot_encoding.confusion.object_alexnet,
             command = encoding.decoding_obj_events.raw_alexnet.only %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               filter(animal_type != "food") %>% 
               mutate(across(c(.obs, .pred), \(x)  x %>% 
                               fct_relevel("dog", "cat", "frog", "spider"))) %>% 
               plot_confusion() +
               this_theme +
               theme(aspect.ratio = 1)
  ),
  tar_target(plot_encoding.decoding_8cat_confusion_flynet.alexnet,
             command = bind_rows(looming = encoding.decoding_8cat_events.endspike_flynet.only,
                                 object = encoding.decoding_8cat_events.endspike_alexnet.only,
                                 .id = "model_type") %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               plot_confusion_8cat(facet_var = model_type) +
               this_theme
  ),
  tar_target(name = plot_encoding.selfreport_combined,
             command = plot_encoding_selfreport_combined(summary_encoding.selfreport %>% 
                                                           filter(data_subset == "animal"),
                                                         summary_encoding.selfreport_pcor %>% 
                                                           filter(data_subset == "animal")) +
               this_theme),
  tar_target(name = plot_encoding.selfreport_flynet.alexnet_pcor,
             command = summary_encoding.selfreport_pcor %>% 
               filter(data_subset == "animal") %>% 
               plot_encoding_selfreport() +
               this_theme),
  tar_target(name = plot_pattern.activation.controlled_flynet,
             command = plot_encoding_controlled_flynet(pattern.activation.controlled_flynet,
                                                       y_var = pattern_expression,
                                                       y_label = "Looming-associated SC activation") +
               this_theme +
               theme(legend.position = "inside",
                     legend.position.inside = c(0,1),
                     legend.justification = c(0,1))),
  tar_target(name = plot_pred.obs.correlation.controlled_flynet,
             command = plot_encoding_controlled_flynet(pattern.activation.controlled_flynet,
                                                       y_var = pred_obs_correlation,
                                                       y_label = "Predicted-actual BOLD correlation") +
               this_theme +
               theme(legend.position = "inside",
                     legend.position.inside = c(0,1),
                     legend.justification = c(0,1))),
  tar_target(name = plot_beta.comparison_object.selfreport,
             command = plot_cormats_encoding_object_selfreport(beta.comparison_category.selfreport) +
               this_theme +
               theme(aspect.ratio = 1)),
  tar_target(name = plot_ggseg.conn_flynet_lat.med,
             command = parcels.wb.connectivity_flynet.only %>% 
               get_parcel_tvals_long() %>% 
               label_parcel_pvals_long(threshold_p = .01) %>% 
               plot_parcel_connectivity_ggseg(fill_col = tval,
                                              pval_col = pval,
                                              ggseg_sides = c("lateral", "medial"),
                                              ggseg_position = hemi ~ side,
                                              max_fill_tval = 10,
                                              viridis_palette = "magma") + 
               theme_brain
  ),
  tar_target(name = plot_ggseg.conn_flynet_dors.vent,
             command = parcels.wb.connectivity_flynet.only %>% 
               get_parcel_tvals_long() %>% 
               label_parcel_pvals_long(threshold_p = .01) %>% 
               plot_parcel_connectivity_ggseg(fill_col = tval,
                                              pval_col = pval,
                                              ggseg_sides = c("dorsal", "ventral"),
                                              max_fill_tval = 10,
                                              viridis_palette = "magma") + 
               theme_brain
  ),
  tar_target(name = plot_ggseg.conn_alexnet_lat.med,
             command = parcels.wb.connectivity_alexnet.only %>% 
               get_parcel_tvals_long() %>% 
               label_parcel_pvals_long(threshold_p = .01) %>% 
               plot_parcel_connectivity_ggseg(fill_col = tval,
                                              pval_col = pval,
                                              ggseg_sides = c("lateral", "medial"),
                                              ggseg_position = hemi ~ side,
                                              max_fill_tval = 10) + 
               theme_brain
  ),
  tar_target(name = plot_ggseg.conn_alexnet_dors.vent,
             command = parcels.wb.connectivity_alexnet.only %>% 
               get_parcel_tvals_long() %>% 
               label_parcel_pvals_long(threshold_p = .01) %>% 
               plot_parcel_connectivity_ggseg(fill_col = tval,
                                              pval_col = pval,
                                              ggseg_sides = c("dorsal", "ventral"),
                                              max_fill_tval = 10) + 
               theme_brain
  ),
  tar_target(name = plot_conn_scatter,
             command =  bind_rows(looming = get_parcel_values_long(parcels.wb.connectivity_flynet.only), 
                                  object = get_parcel_values_long(parcels.wb.connectivity_alexnet.only),
                                  .id = "model_type") %>% 
               relabel_glasser_clt2ggseg() %>% 
               mutate(superregion = case_when(region %in% c("IFJa", "IFJp", "44", "45") ~ "inferior frontal gyrus",
                                              region == "FEF" ~ "frontal eye field",
                                              region %in% c("FFC", "VVC") ~ "ventral occipitotemporal cortex",
                                              region %in% c("IPS1", "MIP", "LIPv", "LIPd", "VIP") ~ "intraparietal cortex",
                                              region == "Amygdala" ~ "amygdala",
                                              region == "Pulv" ~ "pulvinar",
                                              region == "LGN" ~ "LGN",
                                              TRUE ~ NA_character_),
                      superregion = fct_relevel(superregion, 
                                                "amygdala", 
                                                "pulvinar", 
                                                "LGN", 
                                                "ventral occipitotemporal cortex", 
                                                "intraparietal cortex", 
                                                "frontal eye field", 
                                                "inferior frontal gyrus")) %>% 
               filter(!is.na(superregion)) %>% 
               # currently averaging across hemispheres as well
               group_by(model_type, fold_num, superregion) %>% 
               summarize(value = mean(value)) %>% 
               plot_parcel_connectivity_scatter(facet_col = superregion) +
               labs(color = "Connectivity type") +
               this_theme +
               theme(axis.title.y = ggtext::element_markdown(),
                     aspect.ratio = 3)
  ),
  tar_target(name = plot_conn_sig.sim,
             command = wb.connectivity_sig.sim %>% 
               # fine to make the colnames non syntactic, they will get pivoted down into a col before they get sent into plot axes
               rename_with(.fn = \(x) str_replace(x, "_", ": "), .cols = -c(model_type, fold_num)) %>% 
               plot_sig_sim() +
               this_theme +
               theme(aspect.ratio = 3)
  )
)

### plots for supplemental figures, set off for convenience ----
targets_plots_supp <- list(
  tar_target(plot_alexnet.guesses,
             command = get_alexnet_guesses(activations.alexnet_raw,
                                           stims %>% 
                                             select(video, animal_type, looming = has_loom) %>% 
                                             mutate(looming = if_else(looming == 1, "Looming", "No looming")),
                                           imagenet_category_labels,
                                           n_top = 1) %>% 
               plot_predicted_alexnet_categories(n_top = 1) +
               this_theme
  ),
  tar_target(plot_alexnet.guesses.controlled,
             command = get_alexnet_guesses(activations.alexnet_raw_controlled,
                                           metadata_videos_nback %>% 
                                             select(video = filename, animal_type, looming = direction),
                                           imagenet_category_labels,
                                           n_top = 1) %>% 
               plot_predicted_alexnet_categories(n_top = 1) +
               this_theme
  ),
  tar_target(plot_alexnet.guesses.controlled_framewise,
             command = plot_predicted_alexnet_categories_framewise(activations.alexnet_raw_controlled,
                                                                   metadata_videos_nback %>% 
                                                                     select(video = filename, animal_type, looming = direction),
                                                                   imagenet_category_labels,
                                                                   lump_prop = .125)),
  tar_target(plot_auroc.cross_encoding.discrim.looming,
             command = predata_auroc.cross_encoding.discrim.looming %>% 
               roc_curve(truth = animal_value, .pred_outcome_loom) %>% 
               autoplot() +
               labs(title = "classifying object using looming variable",
                    color = "Object category\nbeing classified") +
               this_theme +
               theme(aspect.ratio = 1,
                     legend.position = "inside",
                     legend.position.inside = c(0,1),
                     legend.justification.inside = c(0,1))
             ),
  tar_target(name = plot_auroc.cross_encoding.discrim.object,
             command = predata_auroc.cross_encoding.discrim.object %>% 
               roc_curve(truth = has_loom, .pred_animal) %>% 
               autoplot() +
               labs(title = "classifying looming using object variables",
                    color = "Object category used\nto classify looming") +
               this_theme +
               theme(aspect.ratio = 1,
                     legend.position = "inside",
                     legend.position.inside = c(0,1),
                     legend.justification.inside = c(0,1))
  )
)

### general figure image files ----

targets_figs <- list(
  tar_target(fig_ratings,
             command = ggsave(here::here("ignore", "figs", "naturalistic_ratings.png"),
                              plot = plot_ratings,
                              width = 3000,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_schematic_unit.activation,
             command = ggsave(here::here("ignore", "figs", "naturalistic_schematic_activation.png"),
                              plot = plot_schematic_unit.activation,
                              width = 1500,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_schematic_encoding.bold.pred,
             command = ggsave(here::here("ignore", "figs", "naturalistic_schematic_encoding_bold.png"),
                              plot = plot_schematic_encoding.bold.pred,
                              width = 1500,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_schematic_bold,
             command = ggsave(here::here("ignore", "figs", "naturalistic_schematic_bold.png"),
                              plot = plot_schematic_bold,
                              width = 1500,
                              height = 700,
                              units = "px"),
             format = "file"),
  tar_target(fig_encoding.perf,
             command = ggsave(here::here("ignore", "figs", "naturalistic_encoding_perf.png"),
                              plot = plot_encoding.perf,
                              width = 3500,
                              height = 2000,
                              units = "px"),
             format = "file"),
  tar_target(fig_encoding.perf_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "naturalistic_encoding_perf_flynet.alexnet_sc.png"),
                              plot = plot_encoding.perf_flynet.alexnet +
                                scale_color_manual(values = c("looming" = "green3", "object" = "blue3")),
                              width = 1000,
                              height = 1800,
                              units = "px"),
             format = "file"),
  tar_target(fig_encoding.pcor_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "naturalistic_encoding_pcor_flynet_alexnet.png"),
                              plot = plot_encoding.pcor_flynet.alexnet,
                              width = 3000,
                              height = 2000,
                              units = "px"),
             format = "file")
)

#### SANS 2025 figure image files ----

theme_poster <- theme_bw(base_size = 16) +
  theme(plot.background = element_blank(),
        legend.background = element_blank())

colors_sans2025_loom.noloom <- c("Looming" = "#8ace00", "No looming" = "#527a00")
colors_sans2025_looming.object <- c("looming" = "#8ace00", "object" = "#926372")

targets_figs_sans2025 <- list(
  tar_target(name = fig_sans2025_ratings,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_ratings.png"),
                              plot = plot_ratings_nofood +
                                scale_color_manual(values = colors_sans2025_loom.noloom) +
                                labs(subtitle = NULL,
                                     x = NULL) +
                                theme(legend.position = "inside", 
                                      legend.position.inside = c(1,1), 
                                      legend.justification = c(1,1)),
                              width = 2000,
                              height = 900,
                              units = "px"),
             format = "file"),
  tar_target(fig_sans2025_bold.object_auroc,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_bold_object_auroc.png"),
                              plot = plot_bold.object_auroc +
                                guides(color = guide_legend(override.aes = list(linewidth = 3))) +
                                scale_color_manual(values = colors_sans2025_loom.noloom) +
                                theme(legend.position = "inside",
                                      legend.position.inside = c(1, 0),
                                      legend.justification = c(1, 0)),
                              width = 2000,
                              height = 2300,
                              units = "px"),
             format = "file"),
  tar_target(fig_sans2025_encoding.perf_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_encoding_perf_flynet.alexnet_sc.png"),
                              plot = plot_encoding.perf_flynet.alexnet +
                                scale_color_manual(values = colors_sans2025_looming.object) +
                                # the labels don't need to be angled anymore
                                guides(x = guide_axis()),
                              width = 1800,
                              height = 1600,
                              units = "px"),
             format = "file"),
  tar_target(fig_sans2025_encoding.decoding_confusion_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_encoding_object_confusion_flynet.alexnet_sc.png"),
                              plot = plot_encoding.decoding_confusion_flynet.alexnet +
                                theme(aspect.ratio = 1,
                                      strip.text = element_text(size = rel(1))),
                              width = 4200,
                              height = 2000,
                              units = "px"),
             format = "file"),
  tar_target(fig_sans2025_encoding.selfreport.pcor,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_encoding_selfreport_pcor.png"),
                              plot = plot_encoding.selfreport_flynet.alexnet_pcor +
                                scale_color_manual(values = colors_sans2025_looming.object),
                              width = 1800,
                              height = 800,
                              units = "px"),
             format = "file")
)

#### manuscript 2025 figure image files ----

colors_ms_loom <- "#ff1fbc"
colors_ms_loom.noloom <- c("Looming" = colors_ms_loom, "No looming" = scales::muted(colors_ms_loom, l = 50, c = 30))
colors_ms_looming.object <- c("looming" = colors_ms_loom, "object" = "#28afb0")
colors_ms_looming.object.conj <- c("looming > object" = colors_ms_loom, 
                                   "object > looming" = colors_ms_looming.object[["object"]], # double bracket index to drop old name
                                   "conjunction" = "#9467b6")

targets_figs.ms <- list(
  tar_target(fig_ms_schematic_unit.activation,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_schematic_activation.png"),
                              plot = plot_schematic_unit.activation,
                              width = 1500,
                              height = 1200,
                              units = "px",
                              device = ragg::agg_png),
             format = "file"),
  tar_target(fig_ms_schematic_encoding.bold.pred,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_schematic_encoding_bold.png"),
                              plot = plot_schematic_encoding.bold.pred,
                              width = 1500,
                              height = 1200,
                              units = "px",
                              device = agg_png),
             format = "file"),
  tar_target(fig_ms_schematic_bold,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_schematic_bold.png"),
                              plot = plot_schematic_bold,
                              width = 1500,
                              height = 700,
                              units = "px",
                              device = agg_png),
             format = "file"),
  tar_target(fig_ms_encoding.perf_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_encoding_perf_flynet.alexnet_sc.png"),
                              plot = plot_encoding.perf_flynet.alexnet +
                                scale_color_manual(values = colors_ms_looming.object) +
                                # the labels don't need to be angled anymore
                                guides(x = guide_axis()),
                              width = 1600,
                              height = 1600,
                              units = "px",
                              device = agg_png),
             format = "file"),
  tar_target(fig_ms_encoding.confusion.looming_flynet,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_encoding.confusion.looming.png"),
                              plot = plot_encoding.confusion.looming_flynet,
                              width = 1600,
                              height = 900,
                              units = "px",
                              device = agg_png),
             format = "file"),
  tar_target(fig_ms_encoding.confusion.object_alexnet,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_encoding.confusion.object.png"),
                              plot = plot_encoding.confusion.object_alexnet,
                              width = 2000,
                              height = 1200,
                              units = "px",
                              device = agg_png),
             format = "file"),
  tar_target(fig_ms_pattern.activation.controlled_flynet,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_pattern.activation.controlled_flynet.png"),
                              plot = plot_pattern.activation.controlled_flynet,
                              width = 1600,
                              height = 1400,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_pred.obs.correlation.controlled_flynet,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_pred.obs.correlation.controlled_flynet.png"),
                              plot = plot_pred.obs.correlation.controlled_flynet,
                              width = 1600,
                              height = 1400,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_conn_scatter,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_conn.scatter.png"),
                              plot = plot_conn_scatter + 
                                scale_color_manual(values = colors_ms_looming.object, na.translate = FALSE) +
                                theme(legend.position = "inside",
                                      legend.position.inside = c(1,1),
                                      legend.justification.inside = c(1,1)),
                              width = 2400,
                              height = 2400,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_conn_sig.sim,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_conn.sig.sim.png"),
                              plot = plot_conn_sig.sim + 
                                scale_color_manual(values = colors_ms_looming.object, na.translate = FALSE) +
                                theme(legend.position = "inside",
                                      legend.position.inside = c(0,1),
                                      legend.justification.inside = c(0,1)),
                              width = 2400,
                              height = 2400,
                              units = "px"),
             format = "file"),
  tar_target(name = fig_ms_ratings,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_ratings.png"),
                              plot = plot_ratings_nofood +
                                scale_color_manual(values = colors_ms_loom.noloom) +
                                labs(subtitle = NULL) +
                                theme(legend.position = "inside", 
                                      legend.position.inside = c(0,0), 
                                      legend.justification = c(0,0)),
                              width = 2000,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_encoding.selfreport.pcor,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_encoding_selfreport_pcor.png"),
                              plot = plot_encoding.selfreport_flynet.alexnet_pcor +
                                scale_color_manual(values = colors_ms_looming.object),
                              width = 1800,
                              height = 800,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_encoding.selfreport.combined,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_encoding_selfreport_combined.png"),
                              plot = plot_encoding.selfreport_combined +
                                scale_color_manual(values = colors_ms_looming.object) +
                                scale_fill_manual(values = colors_ms_looming.object),
                              width = 2500,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms_beta.comparison_object.selfreport,
             command = ggsave(here::here("ignore", "figs", "ms_naturalistic_beta.comparison_object.selfreport.png"),
                              plot = plot_beta.comparison_object.selfreport,
                              width = 4200,
                              height = 3600,
                              units = "px"),
             format = "file"),
  tar_target(name = fig_ms_statmaps_pre.mricrogl,
             command = vctrs::vec_c(flynet = statmap.wb.model.connectivity_flynet.only,
                                    alexnet = statmap.wb.model.connectivity_alexnet.only,
                                    conjunction = statmap.wb.model.connectivity_flynet.conj.alexnet,
                                    difference = statmap.wb.model.connectivity_flynet.minus.alexnet,
                                    labels = wb.atlas_selected.rois),
             format = "file")
)

targets_figs.ms.supp <- list(
  tar_target(fig_ms.supp_encoding.confusion.looming,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_encoding.confusion.looming.png"),
                              plot = plot_encoding.confusion.looming,
                              width = 2200,
                              height = 2600,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_encoding.confusion.looming_flynet_events.raw,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_encoding.confusion.looming_events.raw.png"),
                              plot = plot_encoding.confusion.looming_flynet_events.raw,
                              width = 1600,
                              height = 900,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_encoding.confusion.object,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_encoding.confusion.object.png"),
                              plot = plot_encoding.confusion.object,
                              width = 3000,
                              height = 2000,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_auroc.cross_encoding.discrim.looming,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_auroc.cross_encoding.discrim.looming.png"),
                              plot = plot_auroc.cross_encoding.discrim.looming,
                              width = 2000,
                              height = 2000,
                              units = "px")),
  tar_target(fig_ms.supp_auroc.cross_encoding.discrim.object,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_auroc.cross_encoding.discrim.object.png"),
                              plot = plot_auroc.cross_encoding.discrim.object,
                              width = 2000,
                              height = 2000,
                              units = "px")),
  tar_target(fig_ms.supp_alexnet.guesses,
             command = ggsave(here::here("ignore", "figs", "ms.supp_alexnet.guesses.png"),
                              plot = plot_alexnet.guesses +
                                scale_fill_manual(values = colors_ms_loom.noloom) +
                                theme(legend.position = "inside",
                                      legend.position.inside = c(0,1),
                                      legend.justification.inside = c(0,1)),
                              width = 1200,
                              height = 2000,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_alexnet.guesses_controlled,
             command = ggsave(here::here("ignore", "figs", "ms.supp_alexnet.guesses_controlled.png"),
                              plot = plot_alexnet.guesses_controlled +
                                scale_fill_manual(values = colors_ms_loom.noloom) +
                                theme(legend.position = "inside",
                                      legend.position.inside = c(0,1),
                                      legend.justification.inside = c(0,1)),
                              width = 1200,
                              height = 2000,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_norms.controlled,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_norms.controlled.png"),
                              plot = plot_controlled_norms +
                                this_theme +
                                theme(legend.position = "inside", 
                                      legend.position.inside = c(0,0), 
                                      legend.justification = c(0,0)),
                              width = 2000,
                              height = 1200,
                              units = "px"),
             format = "file"),
  tar_target(fig_ms.supp_beta.comparison_object.selfreport,
             command = ggsave(here::here("ignore", "figs", "ms.supp_naturalistic_beta.comparison_object.selfreport.png"),
                              plot = plot_beta.comparison_object.selfreport,
                              width = 4200,
                              height = 2000,
                              units = "px"),
             format = "file")
)

## final combo call ----

# tar_render needs pandoc, which comes with RStudio (Server) but not with any of the R packages
# RStudio Server is only installed on the head node, in a folder that is not mounted on any of the compute nodes
# so when targets runs headless, we need to add to PATH a home-made copy of pandoc in a compute node folder
# rmarkdown::find_pandoc(cache = FALSE, dir = "/home/data/eccolab/Code/renvCache")
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/data/eccolab/Code/renvCache", sep = ":"))

list(
  targets_scripts,
  targets_scripts_naturalistic,
  targets_stimuli,
  targets_stimlists,
  targets_encoding.models,
  targets_qc,
  targets_fmri_across.subject,
  targets_controlled,
  targets_beh,
  targets_plots,
  targets_plots_supp,
  targets_figs,
  targets_figs_sans2025,
  targets_figs.ms,
  targets_figs.ms.supp,
  tar_render(rmd_ms_stats, here::here("code", "R", "fmri_naturalistic_ms_report.rmd"))
)
