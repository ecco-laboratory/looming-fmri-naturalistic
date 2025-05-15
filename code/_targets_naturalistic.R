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
               "qualtRics"), # packages that your targets need to run
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
             "code/R/parse-events.R",
             "code/R/parse-confounds.R",
             "code/R/parse-bold.R",
             "code/R/call-spm.R",
             "code/R/call-canlabtools.R",
             "code/R/make-rdms.R",
             "code/R/utils/tidymodels-helpers.R",
             "code/R/proc-post-canlabtools.R",
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
               arrange(desc(n_spike_runs)))
)

# targets: maps out by subject x run ----

# Use the participants.tsv file as a code-agnostic way of tracking which subjects to use
# It must be edited MANUALLY to label subjects as group "use" once their fmriqc has been checked and approved
# that way, only subjects manually approved will be included in these analyses
participants <- inject(here::here(!!!path_here_fmri, "participants.tsv")) %>% 
  read_tsv() %>% 
  filter(group == "use") %>% 
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

# these take in a single SPM.mat (because )
# the full list of parcels is hard-coded in the matlab_parcellate_betas script
# so changing that script will invalidate this whole target
# the betas are not numbered in a consistent order across subjects
# because they are numbered chronologically by occurrence in the runs
# so it makes sense to have them ALL depend on the level1
# from which we can associate the beta number and the condition name
# and then only map them to files at the end
subtargets_betas.by.parcel <- list(
  # this first target that actually runs the code is NOT format = "file" because the SECOND mapping one is
  # splitting them up instead of using tar_files so that the initial code doesn't re-run EVERY time
  tar_target(name = betas.by.parcel.endspike,
            command = {
              out_path <- file.path(dirname(spm.level1.endspike), "betas_by_parcel")
              dir.create(out_path, showWarnings = FALSE)
              
              matlab_commands = c(
                assign_variable("model_path", spm.level1.endspike),
                assign_variable("out_folder", out_path),
                call_script(matlab_parcellate_betas)
              )
              
              run_matlab_target(matlab_commands, out_path, matlab_path)
            }),
  tar_target(name = betas.by.parcel.map.endspike,
             command = betas.by.parcel.endspike,
             format = "file",
             pattern = map(betas.by.parcel.endspike)),
  tar_target(name = betas.by.parcel.boxcar,
            command = {
              out_path <- file.path(dirname(spm.level1.boxcar), "betas_by_parcel")
              dir.create(out_path, showWarnings = FALSE)
              
              matlab_commands = c(
                assign_variable("model_path", spm.level1.boxcar),
                assign_variable("out_folder", out_path),
                call_script(matlab_parcellate_betas)
              )
              
              run_matlab_target(matlab_commands, out_path, matlab_path)
            }),
  tar_target(name = betas.by.parcel.map.boxcar,
             command = betas.by.parcel.boxcar,
             format = "file",
             pattern = map(betas.by.parcel.boxcar))
)

subtargets_rdms <- list(
  tar_target(name = rdm.endspike,
             command = make_rdms_from_beta(betas.by.parcel.map.endspike),
             pattern = map(betas.by.parcel.map.endspike)),
  tar_target(name = rdm.boxcar,
             command = make_rdms_from_beta(betas.by.parcel.map.boxcar),
             pattern = map(betas.by.parcel.map.boxcar))
)

targets_fmri_by.subject <- make_targets_fmri_by.subject(participants,
                                                        targets_fmri_by.run,
                                                        contrast_names,
                                                        task = "naturalistic",
                                                        additional_targets = c(subtarget_events_combine, 
                                                                               subtargets_encoding.timecourses_combine,
                                                                               subtargets_betas.by.parcel,
                                                                               subtargets_rdms,
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
               unnest(events)),
  tar_target(events.timecourse.boxcartail_all.subs,
             command = events.timecourse_all.subs %>% 
               nest(events = -c(subj_num, fold_num)) %>% 
               mutate(events = map(events, \(x) relabel_timecourse_boxcartail(x))) %>% 
               unnest(events)),
  # 2025-02-18: just for right now make these subcortical only to make them smaller :')
  # 2025-04-03: Phil suggests making RDMs condensed from trialwise to category-wise (either object or object x looming) and looking at classification confusions/clustering them
  tar_combine(name = rdm.endspike_all.subs_subcort,
              targets_fmri_by.subject[["rdm.endspike"]],
              command = list(!!!.x) %>% 
                map(\(x) filter(x, roi == "Bstem_SC")) %>% 
                bind_rows(.id = "target_name") %>% 
                mutate(subj_num = as.integer(str_sub(target_name, start = -4L)),
                       mutate(across(starts_with("condition"), \(x) paste0(str_sub(x, start = -7L), ".mp4")))) %>% 
                filter(condition_row != condition_col) %>% 
                select(-target_name)),
  tar_combine(name = rdm.boxcar_all.subs_subcort,
              targets_fmri_by.subject[["rdm.boxcar"]],
              command = bind_rows(!!!.x, .id = "target_name") %>% 
                filter(roi %in% c("Bstem_SC", "Amygdala")) %>% 
                mutate(subj_num = as.integer(str_sub(target_name, start = -4L))) %>% 
                select(-target_name))
)

### canlabtools-based group analyses so help me god ----

# across-subjects targets, but tar_mapped by ROI
these_rois <- tibble(roi_name = c("amyg", "sc"),
                     roi_canlabtools = c("Amyg", "Bstem_SC"))

subtargets_fmri_canlabtools_by.subject <- tar_map(
  values = participants %>% 
    mutate(bolds = sprintf("bold.smoothed_%s", subject),
           confounds = sprintf("confounds.prespm_%s", subject),
           across(c(bolds, confounds), \(x) syms(str_replace_all(x, "-", "\\.")))),
  # model-based connectivity is set up with single-subject targets to start
  # that depend on the all-subjects pred.encoding.xval target
  # to reduce memory load so that only one subject's whole-brain BOLD needs to be loaded in at a time to calculate the correlations
  tar_target(name = wb.model.connectivity,
             command = canlabtools_fit_model_connectivity(out_path = here::here("ignore", "outputs", sprintf("%s_naturalistic_wb.conn.%s_%s.csv", subject, encoding_type_full, roi_name)),
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
             command = canlabtools_fit_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.%s_%s.csv", encoding_type_full, roi_name)),
                                                    out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.%s_%s.csv", encoding_type_full, roi_name)),
                                                    out_path_betas = here::here("ignore", "outputs", sprintf("naturalistic_betas.%s_%s.csv", encoding_type_full, roi_name)),
                                                    tr_duration = task_defaults_list$tr_duration,
                                                    bolds = bolds_all.subs,
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
  tar_target(name = statmap.perf.encoding,
             command = {
               # canlabtools fmri_data write method appears to forcibly change periods in file names to underscores
               # so we have to feed in something that won't get changed so that the actual output path matches the expected one from here
               this_encoding_label <- str_replace(encoding_type_full, "\\.", "_")
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_perf_%s_%s.nii", this_encoding_label, roi_name)),
                                          roi = roi_canlabtools,
                                          values = summarize_tvals_pre_statmap(perf.encoding.xval),
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_map(
    values = crossing(events_type = c("events.raw", "events.endspike", "events.boxcartail"),
                      n_categories = c("2cat", "5cat", "8cat")) %>% 
      mutate(events_label = case_match(events_type, "events.raw" ~ "", "events.endspike" ~ ".endspike", "events.boxcartail" ~ ".boxcartail"),
             outcome_category = case_match(n_categories, "2cat" ~ "looming", "5cat" ~ "object", "8cat" ~ "object_looming"),
             events_target = syms(sprintf("events.timecourse%s_all.subs", events_label))),
    tar_target(name = encoding.object,
               command = pred.encoding.xval %>%
                 fit_object_by_pattern(events_allsubs = left_join(events_target,
                                                                  beh, 
                                                                  by = c("subj_num", "video_id")),
                                       n_trs_kept_per_run = task_defaults_list$n_trs_kept,
                                       pattern_type = "encoding",
                                       outcome_categories = outcome_category,
                                       n_pls_comp = 20)
    ),
    tar_rep(name = perm.acc_encoding.object,
            command = encoding.object %>% 
              select(preds) %>% 
              unnest(preds) %>% 
              permute_acc_object_by_pattern(n_perms = 100, 
                                            outcome_categories = outcome_category),
            batches = 10,
            reps = 1
    ),
    tar_rep(name = perm.acc_by.object_encoding.object,
            command = encoding.object %>% 
              select(preds) %>% 
              unnest(preds) %>% 
              permute_acc_object_by_pattern(n_perms = 100, 
                                            outcome_categories = outcome_category,
                                            acc_grouping_cols = animal_type),
            batches = 10,
            reps = 1
    ),
    names = c(n_categories, events_type)
  ),
  tar_target(name = encoding.loom_space,
             command = pred.encoding.xval %>% 
               join_beh_to_encoding_trwise(left_join(events.timecourse.endspike_all.subs,
                                                     beh, 
                                                     by = c("subj_num", "video_id"))) %>% 
               select(-starts_with("rating")) %>% 
               # reconstruct run number here for later block permutation by run
               mutate(run_num = (tr_num-1) %/% task_defaults_list$n_trs_kept) %>% 
               # PHIL HAD SET UP HIS PRELIMINARY ANALYSIS TO EXCLUDE FIXATION TIMEPOINTS
               filter(!is.na(animal_type)) %>% 
               fit_model_xval(x_prefix = "voxel",
                              y_prefix = "has_loom",
                              parsnip_model = set_engine(parsnip::pls(mode = "regression", num_comp = 1), engine = "plsr", method = "simpls"),
                              rm_x = TRUE) %>% 
               # the object is too large to save out if we keep all the model fits
               select(preds) %>% 
               unnest(preds)
  ),
  tar_target(name = encoding.selfreport,
             command = pred.encoding.xval %>% 
               get_ratings_by_encoding_space(left_join(events.timecourse.endspike_all.subs,
                                                       beh, 
                                                       by = c("subj_num", "video_id"))) %>% 
               # SO!!! We need to make sure to avoid train-test leakage from the original encoding model to subsequent readout models.
               # we can do this by either:
               # fitting a single overall PLS regression for each of the outcome ratings (I honestly think this is fine. we don't need cross-validated generalization perf)
               # changing this PLS pipeline so that it gets fit along with the main encoding PLS train-test split (I would really prefer not to do this)
               pivot_longer(cols = starts_with("rating"),
                            names_to = "rating_type",
                            values_to = "rating",
                            names_prefix = "rating_") %>% 
               fit_model_2level(nest_vars = c("subj_num", "rating_type"),
                                x_prefix = "voxel",
                                y_var = "rating",
                                # 2025-04-11: Previously, at time of submission for CCN 2025 abstract, had set it to do lm
                                # but I think that is causing the space model to perform way worse on split-half
                                # bc there are more voxels than timepoints. So if you DO want to compare space against time, consider PLS for dimension reduction
                                parsnip_model = parsnip::linear_reg(mode = "regression", engine = "lm"),
                                rm_x = TRUE) %>% 
               mutate(coefs = map(model, \(x) pluck(extract_fit_engine(x), "coefficients"))) %>% 
               select(subj_num, rating_type, coefs, preds)
  ),
  # SEE ABOVE FOR THE DEFINITION OF THE INDIVIDUAL SUBJECT CONNECTIVITY TARGETS
  subtargets_fmri_canlabtools_by.subject,
  tar_combine(name = wb.model.connectivity,
              subtargets_fmri_canlabtools_by.subject[["wb.model.connectivity"]],
              # output as file in case Phil wants to read them in himself
              command = {
                out_path <- here::here("ignore", "outputs", sprintf("naturalistic_wb.conn.%s_%s.csv", encoding_type_full, roi_name))
                
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
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_%s_%s.nii", this_encoding_label, roi_name)),
                                          roi = NULL,
                                          values = tvals.wb.model.connectivity,
                                          threshold_p = .05,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  names = encoding_type_full
)

subtargets_fmri_canlabtools_by.roi <- tar_map(
  values = these_rois %>% 
    mutate(bolds_all.subs = syms(sprintf("bold.masked.%s_all.subs", roi_name))),
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
  tar_map(
    values = crossing(events_type = c("events.raw", "events.endspike", "events.boxcartail"),
                      n_categories = c("2cat", "5cat", "8cat")) %>% 
      mutate(events_label = case_match(events_type, "events.raw" ~ "", "events.endspike" ~ ".endspike", "events.boxcartail" ~ ".boxcartail"),
             outcome_category = case_match(n_categories, "2cat" ~ "looming", "5cat" ~ "object", "8cat" ~ "object_looming"),
             events_target = syms(sprintf("events.timecourse%s_all.subs", events_label))),
    tar_target(name = bold.object,
               command = bolds_all.subs %>%
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
              permute_acc_object_by_pattern(n_perms = 100, 
                                            outcome_categories = outcome_category),
            batches = 10,
            reps = 1
    ),
    tar_target(name = encoding.object_flynet.alexnet,
               command = fit_object_by_pattern(path_pattern_allsubs = pred.encoding.xval_flynet.only_sc,
                                               events_allsubs = left_join(events_target,
                                                                          beh, 
                                                                          by = c("subj_num", "video_id")),
                                               n_trs_kept_per_run = task_defaults_list$n_trs_kept,
                                               path_pattern_allsubs_2 = pred.encoding.xval_alexnet.only_sc,
                                               pattern_type = "encoding",
                                               outcome_categories = outcome_category,
                                               n_pls_comp = 20)
    ),
    tar_rep(name = perm.acc_encoding.object_flynet.alexnet,
            command = encoding.object_flynet.alexnet %>% 
              select(preds) %>% 
              unnest(preds) %>% 
              permute_acc_object_by_pattern(n_perms = 100, 
                                            outcome_categories = outcome_category),
            batches = 10,
            reps = 1
    ),
    names = c(n_categories, events_type)
  ),
  names = roi_name
)

# I don't think these can be tar_mapped more cleanly because they depend on betas and stuff from the base encoding models
subtargets_fmri_canlabtools_compare.models <- tar_map(
  values = these_rois %>% 
    mutate(bolds_all.subs = syms(sprintf("bold.masked.%s_all.subs", roi_name)),
           betas.flynet.onoff = syms(sprintf("betas.encoding.xval_flynet.onoff_%s", roi_name)),
           betas.alexnet.onoff = syms(sprintf("betas.encoding.xval_alexnet.onoff_%s", roi_name)),
           betas.flynet.alexnet = syms(sprintf("betas.encoding.xval_flynet.alexnet_%s", roi_name)),
           preds.flynet = syms(sprintf("pred.encoding.xval_flynet.only_%s", roi_name)),
           preds.alexnet = syms(sprintf("pred.encoding.xval_alexnet.only_%s", roi_name)),
           perf.flynet = syms(sprintf("perf.encoding.xval_flynet.only_%s", roi_name)),
           perf.alexnet = syms(sprintf("perf.encoding.xval_alexnet.only_%s", roi_name)),
           wb.model.connectivity.flynet = syms(sprintf("wb.model.connectivity_flynet.only_%s", roi_name)),
           wb.model.connectivity.alexnet = syms(sprintf("wb.model.connectivity_alexnet.only_%s", roi_name)),
           tvals.wb.model.connectivity.flynet = syms(sprintf("tvals.wb.model.connectivity_flynet.only_%s", roi_name)),
           tvals.wb.model.connectivity.alexnet = syms(sprintf("tvals.wb.model.connectivity_alexnet.only_%s", roi_name))),
  #### performance from joint encoding models with one set of activations knocked out ----
  tar_target(name = encoding.xval_flynet.minus.onoff,
             command = canlabtools_pred_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.flynet.minus.onoff_%s.csv", roi_name)),
                                                     out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.flynet.minus.onoff_%s.csv", roi_name)),
                                                     tr_duration = task_defaults_list$tr_duration,
                                                     bolds = bolds_all.subs,
                                                     activations1 = activations.flynet_all.subs,
                                                     betas = betas.flynet.onoff,
                                                     script = matlab_pred_pls),
             format = "file"),
  tar_target(name = encoding.xval_alexnet.minus.onoff,
             command = canlabtools_pred_encoding_pls(out_path_perf = here::here("ignore", "outputs", sprintf("naturalistic_perf.alexnet.minus.onoff_%s.csv", roi_name)),
                                                     out_path_pred = here::here("ignore", "outputs", sprintf("naturalistic_pred.alexnet.minus.onoff_%s.csv", roi_name)),
                                                     tr_duration = task_defaults_list$tr_duration,
                                                     bolds = bolds_all.subs,
                                                     activations1 = activations.alexnet_all.subs,
                                                     betas = betas.alexnet.onoff,
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
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_perf_flynet_minus_alexnet_%s.nii", roi_name)),
                                          roi = roi_canlabtools,
                                          values = summarize_tvals_pre_statmap(perf.flynet, perf.alexnet),
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_target(name = tvals.wb.model.connectivity_flynet.minus.alexnet,
             command = {
               connectivity.flynet <- vroom::vroom(wb.model.connectivity.flynet, 
                                                   delim = ",", 
                                                   col_names = FALSE, 
                                                   col_types = c(.default = "d"))
               
               connectivity.alexnet <- vroom::vroom(wb.model.connectivity.alexnet, 
                                                    delim = ",", 
                                                    col_names = FALSE, 
                                                    col_types = c(.default = "d"))
               
               summarize_tvals_pre_statmap(in_data_1 = connectivity.flynet,
                                           in_data_2 = connectivity.alexnet,
                                           fun_compare = magrittr::subtract)
             }),
  tar_target(name = tvals.wb.model.connectivity_flynet.conj.alexnet,
             command = {
               # 2025-04-09: Phil wants it this way (Tom Nichols citation?)
               # threshold the individual maps first
               tvals.flynet <- tvals.wb.model.connectivity.flynet %>% 
                 threshold_tvals_pre_statmap(threshold_p = .05)
               tvals.alexnet <- tvals.wb.model.connectivity.alexnet %>% 
                 threshold_tvals_pre_statmap(threshold_p = .05)
               
               pmin(tvals.flynet, tvals.alexnet)

             }),
  tar_target(name = statmap.wb.model.connectivity_flynet.minus.alexnet,
             command = {
               tvals_flynet <- vroom::vroom(wb.model.connectivity.flynet, 
                                                   delim = ",", 
                                                   col_names = FALSE, 
                                                   col_types = c(.default = "d")) %>% 
                 summarize_tvals_pre_statmap() %>% 
                 threshold_tvals_pre_statmap(threshold_p = .05)
               
               tvals_alexnet <- vroom::vroom(wb.model.connectivity.alexnet, 
                                                    delim = ",", 
                                                    col_names = FALSE, 
                                                    col_types = c(.default = "d")) %>% 
                 summarize_tvals_pre_statmap() %>% 
                 threshold_tvals_pre_statmap(threshold_p = .05)
               
                 # there are many voxels where the difference is supra-threshold non-zero but the higher value isn't supra-threshold on its own
                 # so for the diff statmap, zero out voxels ahead of time where the higher value on its own doesn't exceed 0
                 tvals_diff <- tvals.wb.model.connectivity_flynet.minus.alexnet
                 tvals_diff[tvals_diff > 0 & tvals_flynet <= 0] <- 0
                 tvals_diff[tvals_diff < 0 & tvals_alexnet <= 0] <- 0
                 
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_flynet_minus_alexnet_%s.nii", roi_name)),
                                          roi = NULL,
                                          values = tvals_diff,
                                          threshold_p = .05,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  tar_target(name = statmap.wb.model.connectivity_flynet.conj.alexnet,
             command = {
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_flynet_conj_alexnet_%s.nii", roi_name)),
                                          roi = NULL,
                                          values = tvals.wb.model.connectivity_flynet.conj.alexnet,
                                          script = matlab_export_statmap)
             },
             format = "file"),
  names = roi_name
)

#### other downstream-of-canlabtools-encoding targets ----
subtargets_fmri_canlabtools_combined <- list(
  tar_combine(name = perf.encoding_combined,
              tar_select_targets(c(subtargets_fmri_canlabtools_by.roi,
                                   subtargets_fmri_canlabtools_compare.models), 
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
                  separate_wider_delim(cols = target_name, delim = "_", names = c(NA, "encoding_type", "roi"))
              }),
  tar_combine(name = pred.encoding.xval_sc,
              tar_select_targets(c(subtargets_fmri_canlabtools_by.roi,
                                   subtargets_fmri_canlabtools_compare.models), 
                                 starts_with("pred.encoding.xval") & ends_with("sc")),
              command = {
                vctrs::vec_c(!!!.x) %>% 
                  map(\(x) load_encoding_pred_allsubs(x)) %>% 
                  bind_rows(.id = "target_name") %>% 
                  mutate(encoding_type = str_split_i(target_name, "_", -2L)) %>% 
                  select(-target_name)
              }),
  tar_combine(name = pred.encoding.xval_amyg,
              tar_select_targets(c(subtargets_fmri_canlabtools_by.roi,
                                   subtargets_fmri_canlabtools_compare.models), 
                                 starts_with("pred.encoding.xval") & ends_with("amyg")),
              command = {
                vctrs::vec_c(!!!.x) %>% 
                  map(\(x) load_encoding_pred_allsubs(x)) %>% 
                  bind_rows(.id = "target_name") %>% 
                  mutate(encoding_type = str_split_i(target_name, "_", -2L)) %>% 
                  select(-target_name)
              }),
  tar_target(name = pcor.encoding_flynet.alexnet_sc,
             command = calc_pcor_encoding_on_bold(bold.masked.sc_all.subs,
                                                  pred.encoding.xval_flynet.only_sc,
                                                  pred.encoding.xval_alexnet.only_sc,
                                                  "flynet",
                                                  "alexnet")
  ),
  tar_target(name = pcor.encoding_flynet.alexnet_amyg,
             command = calc_pcor_encoding_on_bold(bold.masked.amyg_all.subs,
                                                  pred.encoding.xval_flynet.only_amyg,
                                                  pred.encoding.xval_alexnet.only_amyg,
                                                  "flynet",
                                                  "alexnet")
  )
)

## targets imported from the controlled task ----
# these do NOT auto update when new subjects are processed for naturalistic
# you must tar_make the original controlled targets for that tracking
#  but if you do tar_make controlled, then the correct naturalistic dependency targets will be marked as out of date

metadata_videos_nback <- tar_read(metadata_videos_nback, store = here::here("ignore", "_targets", "controlled"))

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
  tar_target(pattern.expression.controlled_flynet_sc,
             command = calc_controlled_pattern_expression(cons_level1.5_controlled,
                                                          activations.flynet_raw_controlled,
                                                          metadata_videos_nback,
                                                          betas.encoding.xval_flynet.only_sc)),
  tar_target(pattern.expression.controlled_alexnet_sc,
             command = calc_controlled_pattern_expression(cons_level1.5_controlled,
                                                          activations.alexnet_raw_controlled,
                                                          metadata_videos_nback,
                                                          betas.encoding.xval_alexnet.only_sc))
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
                 # grand mean center it
                 mutate(animal_type = factor(animal_type, levels = c("dog", "cat", "frog", "spider"))) %>% 
                 rename(this_rating = paste0("rating_", rating_type)) %>% 
                 group_by(subj_num, animal_type, has_loom) %>% 
                 summarize(rating_mean = mean(this_rating)) %>% 
                 lm(rating_mean ~ has_loom * animal_type, data = .) %>% 
                 anova())
  ),
  tar_target(name = summary_perf.encoding,
             command = perf.encoding_combined %>% 
               group_by(roi, encoding_type) %>% 
               mutate(r2 = perf^2) %>% 
               summarize(r_mean = mean(perf),
                         r_se = sd(perf)/sqrt(length(perf)),
                         r2_mean = mean(r2),
                         r2_se = sd(r2)/sqrt(length(r2)),
                         .groups = "drop") %>% 
               mutate(r_ci95_lower = r_mean - 1.96*r_se,
                      r_ci95_upper = r_mean + 1.96*r_se,
                      r2_ci95_lower = r2_mean - 1.96*r2_se,
                      r2_ci95_upper = r2_mean + 1.96*r2_se)
  ),
  tar_target(name = summary_perf.encoding_alexnet.minus.flynet,
             command = perf.encoding_combined %>% 
               filter(encoding_type %in% c("flynet.only", "alexnet.only")) %>% 
               pivot_wider(names_from = encoding_type, values_from = perf) %>% 
               # a priori we know alexnet is higher so this will get it signed positive. whatever
               mutate(perf_diff = alexnet.only - flynet.only,
                      r2_diff = alexnet.only^2 - flynet.only^2) %>% 
               group_by(roi) %>% 
               summarize(diff_mean = mean(perf_diff),
                         diff_se = sd(perf_diff)/sqrt(length(perf_diff)),
                         diff_cohens.d = (mean(alexnet.only) - mean(flynet.only)) / sqrt((var(alexnet.only) + var(flynet.only))/2),
                         r2.diff_mean = mean(r2_diff),
                         r2.diff_cohens.d = (mean(alexnet.only^2) - mean(flynet.only^2)) / sqrt((var(alexnet.only^2) + var(flynet.only^2))/2))
  ),
  tar_target(name = summary_pcor.encoding_flynet.alexnet_sc,
             command = pcor.encoding_flynet.alexnet_sc %>% 
               group_by(term) %>% 
               summarize(pcor_mean = mean(estimate),
                         pcor_se = sd(estimate)/sqrt(length(estimate))) %>% 
               mutate(pcor_ci95_lower = pcor_mean - 1.96*pcor_se,
                      pcor_ci95_upper = pcor_mean + 1.96*pcor_se)
             ),
  tar_target(name = encoding.selfreport_flynet.alexnet.cbound_sc,
             command = get_ratings_by_encoding_space(path_pred_allsubs = pred.encoding.xval_flynet.only_sc,
                                                     events_allsubs = left_join(events.timecourse.endspike_all.subs,
                                                                                beh, 
                                                                                by = c("subj_num", "video_id")),
                                                     path_pred_allsubs_2 = pred.encoding.xval_alexnet.only_sc) %>% 
               # SO!!! We need to make sure to avoid train-test leakage from the original encoding model to subsequent readout models.
               # we can do this by either:
               # fitting a single overall PLS regression for each of the outcome ratings (I honestly think this is fine. we don't need cross-validated generalization perf)
               # changing this PLS pipeline so that it gets fit along with the main encoding PLS train-test split (I would really prefer not to do this)
               pivot_longer(cols = starts_with("rating"),
                            names_to = "rating_type",
                            values_to = "rating",
                            names_prefix = "rating_") %>% 
               fit_model_2level(nest_vars = c("subj_num", "rating_type"),
                                x_prefix = "voxel",
                                y_var = "rating",
                                # 2025-04-11: Previously, at time of submission for CCN 2025 abstract, had set it to do lm
                                # but I think that is causing the space model to perform way worse on split-half
                                # bc there are more voxels than timepoints. So if you DO want to compare space against time, consider PLS for dimension reduction
                                parsnip_model = parsnip::linear_reg(mode = "regression", engine = "lm"),
                                rm_x = TRUE) %>% 
               mutate(coefs = map(model, \(x) pluck(extract_fit_engine(x), "coefficients"))) %>% 
               select(subj_num, rating_type, coefs, preds)
  ),
  tar_target(name = encoding.selfreport_rbound_sc,
             command = bind_rows(alexnet = encoding.selfreport_alexnet.only_sc,
                                 flynet = encoding.selfreport_flynet.only_sc,
                                 cbound = encoding.selfreport_flynet.alexnet.cbound_sc,
                                 .id = "model_type") %>% 
               select(subj_num, rating_type, model_type, preds) %>% 
               unnest(preds) %>% 
               # Clearly some code generalization failure thing is making this column of all NAs pop up. idk man
               select(-rating)),
  tar_target(name = summary_encoding.selfreport_sc,
             command = encoding.selfreport_rbound_sc %>% 
               nest(.by = c(model_type, rating_type)) %>% 
               mutate(data = map(data, \(x) group_bootstraps(x, group = subj_num, times = 1000))) %>%
               unnest(data) %>% 
               mutate(cor.overall = map_dbl(splits, \(x) x %>% 
                                              analysis() %$%
                                              cor(.obs, .pred, method = "spearman"),
                                            .progress = "bootstrapping pred-obs correlations"),
                      cor.food = map_dbl(splits, \(x) x %>% 
                                           analysis() %>% 
                                           filter(animal_type == "food") %$%
                                           cor(.obs, .pred, method = "spearman"),
                                         .progress = "bootstrapping pred-obs correlations FOOD ONLY")) %>%
               group_by(rating_type, model_type) %>% 
               summarize(across(starts_with("cor"), 
                                summary_funs_bootstrap))
  ),
  tar_target(name = summary_encoding.selfreport_sc_pcor,
             command = encoding.selfreport_flynet.alexnet_sc %>% 
               pivot_wider(names_from = model_type, values_from = .pred) %>% 
               nest(.by = rating_type) %>% 
               mutate(data = map(data, \(x) group_bootstraps(x, group = subj_num, times = 1000))) %>%
               unnest(data) %>% 
               mutate(coefs.overall = map(splits, \(x) x %>% 
                                            analysis() %>% 
                                            # 2025-04-05: SPEARMAN CORRELATION because the split-half models are (predictably) way noisier out of sample. so there are some wild outliers
                                            mutate(alexnet = dense_rank(alexnet),
                                                   flynet = dense_rank(flynet)) %>% 
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
               select(-splits) %>% 
               pivot_longer(cols = starts_with("coefs"), names_to = "data_subset", values_to = "coefs", names_prefix = "coefs.") %>% 
               unnest(coefs) %>% 
               filter(term != "(Intercept)") %>% 
               group_by(rating_type, data_subset, term) %>% 
               rename(pcor = estimate) %>% 
               summarize(across(pcor, summary_funs_bootstrap),
                         .groups = "drop")
  ),
  tar_target(name = summary_encoding.discrim.looming_sc,
             command = {
               flynet <- calc_perm_pval_object_by_pattern(preds = encoding.object_2cat_events.endspike_flynet.only_sc, 
                                                          perms = perm.acc_encoding.object_2cat_events.endspike_flynet.only_sc)
               alexnet <- calc_perm_pval_object_by_pattern(preds = encoding.object_2cat_events.endspike_alexnet.only_sc, 
                                                           perms = perm.acc_encoding.object_2cat_events.endspike_alexnet.only_sc)
               bind_rows(flynet = flynet,
                         alexnet = alexnet,
                         .id = "encoding_type")
               }),
  tar_target(name = summary_encoding.discrim.looming_by.category_sc,
             command = {
               flynet <- calc_perm_pval_object_by_pattern(preds = encoding.object_2cat_events.endspike_flynet.only_sc, 
                                                          perms = perm.acc_by.object_encoding.object_2cat_events.endspike_flynet.only_sc,
                                                          grouping_cols = animal_type)
               alexnet <- calc_perm_pval_object_by_pattern(preds = encoding.object_2cat_events.endspike_alexnet.only_sc, 
                                                           perms = perm.acc_by.object_encoding.object_2cat_events.endspike_alexnet.only_sc,
                                                           grouping_cols = animal_type)
               bind_rows(flynet = flynet,
                         alexnet = alexnet,
                         .id = "encoding_type")
             }),
  tar_target(name = summary_encoding.discrim.object_sc,
             command = {
               flynet <- calc_perm_pval_object_by_pattern(preds = encoding.object_5cat_events.endspike_flynet.only_sc, 
                                                          perms = perm.acc_encoding.object_5cat_events.endspike_flynet.only_sc)
               alexnet <- calc_perm_pval_object_by_pattern(preds = encoding.object_5cat_events.endspike_alexnet.only_sc, 
                                                           perms = perm.acc_encoding.object_5cat_events.endspike_alexnet.only_sc)
               bind_rows(flynet = flynet,
                         alexnet = alexnet,
                         .id = "encoding_type")
             })
)

targets_fmri_across.subject <- make_targets_fmri_across.subject(targets_fmri_by.subject,
                                                                contrast_names,
                                                                task = "naturalistic",
                                                                additional_targets = c(subtargets_fmri_across.subject,
                                                                                       subtargets_fmri_canlabtools_by.roi,
                                                                                       subtargets_fmri_canlabtools_compare.models,
                                                                                       subtargets_fmri_canlabtools_combined,
                                                                                       subtargets_fmri_summary))

## targets: splitting by parcel ROI (for whole-brainy analyses) ----

if (FALSE) {
  targets_wholebrain <- list(
    tar_eval(
      tar_target(
        name = target_name,
        command = {
          list(!!!input_names) %>% 
            set_names(combine_vals) %>% 
            bind_rows(.id = "subj_num") %>% 
            mutate(subj_num = as.integer(str_sub(subj_num, start = -4L))) %>% 
            separate_wider_delim(cols = condition_row, 
                                 delim = "_", 
                                 names = c("animal1", "loom1", "video11", "video21")) %>% 
            separate_wider_delim(cols = condition_col, 
                                 delim = "_", 
                                 names = c("animal2", "loom2", "video12", "video22")) %>% 
            unite("video1", video11, video21) %>% 
            unite("video2", video12, video22) %>% 
            filter(video1 != video2) %>% 
            mutate(video_sort = map2_chr(video1, video2, \(x, y) paste(sort(c(x, y)), collapse = " "))) %>% 
            distinct(subj_num, roi, video_sort, .keep_all = TRUE) %>% 
            select(-video_sort) %>% 
            mutate(across(starts_with("loom"), \(x) as.integer(str_sub(x, start = -1L))),
                   same_animal = as.integer(animal1 == animal2), 
                   same_loom = as.integer(loom1 == loom2),
                   # per Gower 1966, Gower & Legendre 1996, Solo 2019
                   # Pearson "distance" in this way has a range of 2 if negative correlations are considered far
                   # sqrt(1 - r) is metric, in that it satisfies the triangle inequality for any 3 observations
                   distance = sqrt(1 - correlation))
        }
      ),
      values = map_values_across.task %>% 
        filter(endsWith(suffix, "naturalistic")) %>% 
        distinct(suffix, model_type, combine_vals) %>% 
        mutate(input_names = pmap(list(model_type, suffix, combine_vals), \(a, b, c) syms(paste("rdms_smoothed.4mm", a, b, c, sep = "_"))),
               target_name = syms(paste("rdms_smoothed.4mm", model_type, suffix, sep = "_"))) %>% 
        select(target_name, input_names, combine_vals)
    ),
    tar_target(
      rdm_preplot_endspike_subcort_naturalistic,
      command = {
        rdms_smoothed.4mm_endspike_task.naturalistic %>% 
          filter(roi %in% c("Bstem_SC", "Amygdala")) %>% 
          mutate(across(starts_with("video"), \(x) paste0(x, ".mp4"))) %>% 
          left_join(rdms_beh, by = c("subj_num", "video1", "video2")) %>% 
          # the sorting wasn't the same I guess??
          left_join(rdms_beh, by = c("subj_num", "video1" = "video2", "video2" = "video1")) %>% 
          mutate(diff_pleasantness = coalesce(diff_pleasantness.x, diff_pleasantness.y), 
                 diff_arousal = coalesce(diff_arousal.x, diff_arousal.y), 
                 diff_fear = coalesce(diff_fear.x, diff_fear.y),
                 # so the betas are in units of one slider step
                 across(starts_with("diff"), \(x) x / 10)) %>% 
          select(-ends_with(".x"), -ends_with(".y")) %>% 
          group_by(roi, video1, animal1, loom1, video2, animal2, loom2) %>% 
          summarize(correlation = mean(correlation), 
                    distance = mean(distance), 
                    same_animal = mean(same_animal),
                    same_loom = mean(same_loom),
                    .groups = "drop") %>% 
          mutate(intxn = same_animal * same_loom)
      }
    )
  )
}

## targets: behavioral data ----

targets_beh <- list(
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

targets_beh_old <- list(
  # from 3 colleagues I was able to shake down for ratings in March 2024, lol
  tar_files(name = norms_inperson_raw,
            command = list.files(here::here("ignore", "data", "norm"),
                                 recursive = TRUE,
                                 full.names = TRUE) %>% 
              .[!grepl("debug", .)]
  ),
  tar_target(name = norms_inperson,
             command = norms_inperson_raw %>% 
               read_csv(na = c("", "NA", "None")),
             pattern = map(norms_inperson_raw)),
  tar_target(
    name = norms_qualtrics,
    command = {
      annotations_naturalistic <- read_csv(annotations_videos)
      annotations_controlled <- metadata_videos_nback %>% 
        left_join(read_csv(qualtrics.ids_videos_controlled), by = "filename") %>% 
        # to get the cols in the same format as the naturalistic annotations
        mutate(has_loom = as.numeric(direction == "looming")) %>% 
        select(video_id = filename, qualtrics_id, animal = animal_type, has_loom, hemifield)
      
      annotations <- bind_rows(naturalistic = annotations_naturalistic, 
                               controlled = annotations_controlled,
                               .id = "stim_type")
        
      
      # TODO: Make sure your qualtrics API token is reproducible yet safe from prying eyes
      # TODO also: set up another target so that targets will think the survey is updated when it is
      qualtrics_data <- fetch_survey("SV_6nZOkqUU0MEAnMa",
                                     include_display_order = FALSE)
      
      video_qualtrics_ids <- qualtrics_data %>% 
        select(contains("caption")) %>% 
        sjlabelled::get_label() 
      
      df_video_qualtrics_ids <- tibble(video_num = names(video_qualtrics_ids),
                                       qualtrics_id = video_qualtrics_ids %>% 
                                         str_split_i(pattern = " ", i = 1)) %>% 
        separate_wider_regex(cols = video_num,
                             patterns = c(video_num = "[[:digit:]]+",
                                          "_[[:alpha:]]+",
                                          block_num = "[[:digit:]]+")) %>% 
        mutate(across(ends_with("num"), as.integer))
      
      qualtrics_data %>% 
        filter(Status != "Survey Preview", participantId_check != "test") %>% 
        select(-contains("quit"), -starts_with("Recipient"), -DistributionChannel, -UserLanguage) %>% 
        pivot_longer(cols = c(contains("caption"), 
                              contains("pleasantness"), 
                              contains("arousal"), 
                              contains("fear")), 
                     names_to = c("video_num", ".value", "block_num"), 
                     names_pattern = "([[:digit:]]+)_([[:alpha:]]+)([[:digit:]]+)",
                     names_transform = list(video_num = as.integer,
                                            block_num = as.integer)) %>% 
        left_join(df_video_qualtrics_ids, by = c("video_num", "block_num")) %>% 
        left_join(annotations,
                  by = "qualtrics_id")
    }
  )
)

## targets: plots for showing ----

this_theme <- theme_bw(base_size = 16) +
  theme(plot.background = element_blank(),
        legend.background = element_blank())

targets_plots <- list(
  tar_target(plot_ratings,
             command = plot_selfreport_ratings(beh) +
                 this_theme
  ),
  tar_target(plot_ratings_nofood,
             command = beh %>% 
               filter(animal_type != "food") %>% 
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
             command = plot_sample_timecourse_encoding(pred.encoding.xval_sc,
                                                       encoding_models = c("looming" = "flynet.only",
                                                                           "object" = "alexnet.only")) +
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
             command = plot_auroc_8cat(bold.object_8cat_events.endspike_sc) +
               this_theme
  ),
  tar_target(plot_encoding.perf,
             command = plot_encoding_performance(perf.encoding_combined) +
               this_theme
  ),
  tar_target(plot_encoding.perf_flynet.alexnet_sc,
             command = plot_encoding_performance(perf.encoding_combined, 
                                                 encoding_types = c("looming" = "flynet.only", 
                                                                    "object" = "alexnet.only"), 
                                                 rois = "sc") +
               guides(color = "none") + 
               labs(x = NULL,
                    subtitle = NULL) +
               this_theme +
               theme(axis.title = element_text(size = rel(1.3)),
                     axis.text = element_text(size = rel(1.05)))
  ),
  tar_target(plot_encoding.pcor_flynet.alexnet,
             command = bind_rows(sc = pcor.encoding_flynet.alexnet_sc,
                                 amyg = pcor.encoding_flynet.alexnet_amyg,
                                 .id = "roi") %>% 
               plot_encoding_performance_pcor() +
               this_theme
  ),
  tar_target(plot_encoding.object_confusion_flynet.alexnet,
             command = bind_rows(looming = encoding.object_8cat_events.endspike_flynet.only_sc,
                                 object = encoding.object_8cat_events.endspike_alexnet.only_sc,
                                 .id = "model_type") %>% 
               select(-coefs) %>% 
               unnest(preds) %>% 
               plot_confusion_8cat(facet_var = model_type) +
               this_theme
  ),
  tar_target(name = plot_encoding.selfreport_flynet.alexnet_pcor,
             command = plot_encoding_selfreport_pcor(summary_encoding.selfreport_sc_pcor) +
               this_theme)
)

targets_plots.rdm <- list(
  tar_target(
    plot_rdm_subcort,
    command = {
      preplot <- rdm_preplot_endspike_subcort %>% 
        # stupid shit to re-complete both triangle halves of the square
        bind_rows(rdm_preplot_endspike_subcort %>% 
                    rename(video3 = video1, animal3 = animal1, loom3 = loom1) %>% 
                    rename(video1 = video2, animal1 = animal2, loom1 = loom2) %>% 
                    rename(video2 = video3, animal2 = animal3, loom2 = loom3))
      
      square_bounds_base <- stims %>% 
        # FOR SANS 2025 ABSTRACT WE AREN'T REPORTING FOOD
        filter(animal_type != "food") %>% 
        mutate(animal_type = fct_relevel(animal_type, "dog", "cat", "frog")) %>% 
        arrange(animal_type, has_loom)
      
      square_bounds_intxn <- square_bounds_base %>% 
        count(animal_type, has_loom) %>% 
        mutate(xmax = cumsum(n) + 0.5, 
               xmin = coalesce(lag(xmax), 0.5),
               label = if_else(has_loom == 0,
                               paste("non-looming", animal_type),
                               paste("looming", animal_type)))
      
      preplot %>% 
        mutate(across(starts_with("animal"), \(x) as.integer(factor(x, levels = c("dog", "cat", "frog", "spider", "food")))),
               roi = case_match(roi,
                                "Amygdala" ~ "amygdala",
                                "Bstem_SC" ~ "superior colliculus")) %>% 
        unite(col = "intxn1", animal1, loom1, remove = FALSE) %>% 
        unite(col = "intxn2", animal2, loom2, remove = FALSE) %>% 
        # FOR SANS 2025 ABSTRACT WE AREN'T REPORTING FOOD
        filter(animal1 != 5, animal2 != 5) %>% 
        # set the aesthetics in each layer bc the data don't all have the same cols
        ggplot() + 
        geom_raster(aes(x = fct_reorder(video1, intxn1, .fun = \(x) sort(unique(x))), 
                        y = fct_reorder(video2, intxn2, .fun = \(x) sort(unique(x))),
                        fill = distance)) + 
        geom_rect(aes(xmin = xmin,
                      xmax = xmax,
                      ymin = xmin,
                      ymax = xmax),
                  data = square_bounds_intxn,
                  color = "white",
                  alpha = 0) +
        geom_label(aes(x = xmin, y = xmax, label = label),
                   data = square_bounds_intxn,
                   hjust = 0, vjust = 1,
                   size = 3) +
        facet_wrap(~roi) + 
        scale_fill_viridis_c() + 
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 16) +
        theme(plot.background = element_blank(),
              legend.background = element_blank(),
              aspect.ratio = 1,
              axis.text = element_blank(),
              axis.ticks = element_blank())
    }
  ),
  tar_target(
    plot_spaghetti.intxn,
    command = rdms_smoothed.4mm_endspike_task.naturalistic %>% 
      filter(roi %in% c("Bstem_SC", "Amygdala")) %>% 
      mutate(intxn = if_else(same_animal * same_loom == 1, 
                             "Same object & loom\n(integration)", 
                             "All other pairs"), 
             roi = case_match(roi,
                              "Amygdala" ~ "amygdala",
                              "Bstem_SC" ~ "superior colliculus")) %>% 
      group_by(roi, subj_num, intxn) %>% 
      summarize(distance = mean(distance), .groups = "drop") %>% 
      ggplot(aes(x = intxn, y = distance)) + 
      geom_line(aes(group = subj_num)) + 
      geom_point() + 
      guides(x = guide_axis(angle = 20)) + 
      facet_wrap(~roi) + 
      labs(x = "Pair type", y = "fMRI pattern distance") + 
      theme_bw(base_size = 16) + 
      theme(plot.background = element_blank())
  )
)

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
  tar_target(fig_encoding.perf_flynet.alexnet_sc,
             command = ggsave(here::here("ignore", "figs", "naturalistic_encoding_perf_flynet.alexnet_sc.png"),
                              plot = plot_encoding.perf_flynet.alexnet_sc +
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
  tar_target(fig_sans2025_encoding.perf_flynet.alexnet_sc,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_encoding_perf_flynet.alexnet_sc.png"),
                              plot = plot_encoding.perf_flynet.alexnet_sc +
                                scale_color_manual(values = colors_sans2025_looming.object) +
                                # the labels don't need to be angled anymore
                                guides(x = guide_axis()),
                              width = 1800,
                              height = 1600,
                              units = "px"),
             format = "file"),
  tar_target(fig_sans2025_encoding.object_confusion_flynet.alexnet,
             command = ggsave(here::here("ignore", "figs", "sans2025_naturalistic_encoding_object_confusion_flynet.alexnet_sc.png"),
                              plot = plot_encoding.object_confusion_flynet.alexnet +
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


targets_figs.rdm <- list(
  tar_target(
    fig_rdm_subcort,
    command = ggsave(here::here("ignore", "figs", "rdm_subcort_naturalistic.png"),
                     plot = plot_rdm_subcort,
                     width = 3000,
                     height = 1600,
                     units = "px")
  ),
  tar_target(
    fig_spaghetti.intxn,
    command = ggsave(here::here("ignore", "figs", "spaghetti.intxn_naturalistic.png"),
                     plot = plot_spaghetti.intxn,
                     width = 1600,
                     height = 1600,
                     units = "px")
  )
)

## final combo call ----

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
  targets_figs,
  targets_figs_sans2025
)
