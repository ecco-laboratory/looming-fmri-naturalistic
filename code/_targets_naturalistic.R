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
               "tidyverse",
               "magrittr",
               "glue",
               "rlang",
               "qualtRics"), # packages that your targets need to run
  controller = crew.cluster::crew_controller_slurm(
    workers = 12,
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

## targets: looming stimuli of various kinds ----

folder_videos <- here::here("ignore", "stimuli", "videos", "naturalistic")

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

# targets: maps out by subject x run ----

# Use the participants.tsv file as a code-agnostic way of tracking which subjects to use
# It must be edited MANUALLY to label subjects as group "use" once their fmriqc has been checked and approved
# that way, only subjects manually approved will be included in these analyses
participants <- inject(here::here(!!!path_here_fmri, "participants.tsv")) %>% 
  read_tsv() %>% 
  filter(group == "use") %>% 
  select(subject = participant_id)


subtargets_encoding.timecourses_by.run <- list(
  tar_target(events.timecourse,
             command = make_condition_timecourse(onsets = events,
                                                 tr_duration = task_defaults_list$tr_duration,
                                                 n_trs = task_defaults_list$n_trs_kept)),
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
                    "frog",
                    "spider",
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
  tar_combine(name = events.timecourse,
              targets_fmri_by.run[["events.timecourse"]],
              command = {
                bind_rows(!!!.x, .id = "target_name") %>% 
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
  # 2025-02-18: just for right now make these subcortical only to make them smaller :')
  tar_combine(name = rdm.endspike_all.subs_subcort,
              targets_fmri_by.subject[["rdm.endspike"]],
              command = bind_rows(!!!.x, .id = "target_name") %>% 
                filter(roi %in% c("Bstem_SC", "Amygdala")) %>% 
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
  tar_target(name = encoding.selfreport.xval_time,
             command = pred.encoding.xval %>% 
               get_ratings_by_encoding_time(left_join(events.timecourse.endspike_all.subs,
                                                      beh, 
                                                      by = c("subj_num", "video_id"))) %>% 
               fit_pls_xval(x_prefix = "tr",
                            num_comp = 3)
  ),
  tar_target(name = encoding.selfreport.xval_space,
             command = pred.encoding.xval %>% 
               get_ratings_by_encoding_space(left_join(events.timecourse.endspike_all.subs,
                                                       beh, 
                                                       by = c("subj_num", "video_id"))) %>% 
               fit_pls_xval(x_prefix = "tr",
                            num_comp = 3)
  ),
  # this one is THIRSTY FOR MEMORY bc it requires all subjects' whole brain full timecourses. 172 GB seems to work with 29 subjects
  tar_target(name = wb.model.connectivity,
             command = canlabtools_fit_model_connectivity(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb.conn.%s_%s.csv", encoding_type_full, roi_name)),
                                                          tr_duration = task_defaults_list$tr_duration,
                                                          trs_to_use = 1:task_defaults_list$n_trs_kept + (task_defaults_list$disdaq_duration %/% task_defaults_list$tr_duration),
                                                          bolds = bold.smoothed_all.subs,
                                                          confounds = confounds.prespm_all.subs,
                                                          pred.encoding.roi = pred.encoding.xval,
                                                          script = matlab_fit_model_connectivity),
             format = "file"),
  tar_target(name = statmap.wb.model.connectivity,
             command = {
               connectivity <- vroom::vroom(wb.model.connectivity, delim = ",", col_names = FALSE, col_types = c(.default = "d"))
               # canlabtools fmri_data write method appears to forcibly change periods in file names to underscores
               # so we have to feed in something that won't get changed so that the actual output path matches the expected one from here
               this_encoding_label <- str_replace(encoding_type_full, "\\.", "_")
               canlabtools_export_statmap(out_path = here::here("ignore", "outputs", sprintf("naturalistic_wb_conn_%s_%s.nii", this_encoding_label, roi_name)),
                                          roi = NULL,
                                          values = summarize_tvals_pre_statmap(perf.encoding.xval),
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
           preds.alexnet = syms(sprintf("pred.encoding.xval_alexnet.only_%s", roi_name))),
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

targets_fmri_across.subject <- make_targets_fmri_across.subject(targets_fmri_by.subject,
                                                                contrast_names,
                                                                task = "naturalistic",
                                                                additional_targets = c(subtargets_fmri_across.subject,
                                                                                       subtargets_fmri_canlabtools_by.roi,
                                                                                       subtargets_fmri_canlabtools_compare.models,
                                                                                       subtargets_fmri_canlabtools_combined))

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
                              pleasantness1 = pleasantness_rating, 
                              arousal1 = arousal_rating, 
                              fear1 = fear_rating), 
                      nesting(video2 = video_id, 
                              pleasantness2 = pleasantness_rating, 
                              arousal2 = arousal_rating, 
                              fear2 = fear_rating)) %>% 
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
  tar_target(plot_schematic_encoding.bold.pred,
             command = plot_sample_timecourse_encoding(pred.encoding.xval_sc,
                                                       encoding_models = c("flynet.only",
                                                                           "alexnet.only",
                                                                           "onoff.only")) +
               this_theme +
               theme(aspect.ratio = 1/3)
  ),
  tar_target(plot_schematic_bold,
             command = plot_sample_timecourse_bold(bold.masked.sc_sub.0001) +
               this_theme +
               theme(aspect.ratio = 1/3)
  ),
  tar_target(plot_encoding.perf,
             command = plot_encoding_performance(perf.encoding_combined) +
               this_theme
  ),
  tar_target(plot_encoding.pcor_flynet.alexnet,
             command = bind_rows(sc = pcor.encoding_flynet.alexnet_sc,
                                 amyg = pcor.encoding_flynet.alexnet_amyg,
                                 .id = "roi") %>% 
               plot_encoding_performance_pcor() +
               this_theme
  )
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
  tar_target(fig_schematic_encoding.bold.pred,
             command = ggsave(here::here("ignore", "figs", "naturalistic_schematic_encoding_bold.png"),
                              plot = plot_schematic_encoding.bold.pred,
                              width = 2800,
                              height = 2400,
                              units = "px"),
             format = "file"),
  tar_target(fig_schematic_bold,
             command = ggsave(here::here("ignore", "figs", "naturalistic_schematic_bold.png"),
                              plot = plot_schematic_bold,
                              width = 1500,
                              height = 800,
                              units = "px"),
             format = "file"),
  tar_target(fig_encoding.perf,
             command = ggsave(here::here("ignore", "figs", "naturalistic_encoding_perf.png"),
                              plot = plot_encoding.perf,
                              width = 3500,
                              height = 2000,
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
  targets_stimuli,
  targets_stimlists,
  targets_encoding.models,
  targets_fmri_across.subject,
  targets_beh,
  targets_plots,
  targets_figs
)
