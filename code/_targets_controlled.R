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
  packages = c("tarchetypes",
               "withr",
               "matlabr",
               "tidyverse",
               "magrittr",
               "glue",
               "rlang",
               "qualtRics"), # packages that your targets need to run
  controller = crew.cluster::crew_controller_slurm(
    workers = 8,
    seconds_idle = 15,
    options_cluster = crew.cluster::crew_options_slurm(
      verbose = TRUE,
      script_lines = "#SBATCH --account=default",
      log_output = "/home/%u/log/crew_log_%A.out",
      log_error = "/home/%u/log/crew_log_%A.err",
      # single subject level 1s require 32 GB to run more than snail's pace?
      memory_gigabytes_required = 16,
      cpus_per_task = 1,
      time_minutes = 1339,
      partition = "day-long"
    )
  )
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source(c("code/R/utils/",
             "code/R/parse-demos.R",
             "code/R/make-stimlist.R",
             "code/R/call-activations.R",
             "code/R/define-targets-fmri.R",
             "code/R/parse-qualtrics.R",
             "code/R/parse-events.R",
             "code/R/parse-confounds.R",
             "code/R/parse-bold.R",
             "code/R/call-spm.R",
             "code/R/call-canlabtools.R",
             "code/R/make-plots.R"))
# Regular source this script because it's not called by a target per se
# but is necessary for target construction
# source("code/R/set-study-defaults.R")

## other useful global vars ----

# fMRI parameters that pass into multiple targets
nback_rating_duration <- 12 # tr units
nback_blocks_per_run <- 4
nback_block_threatenings <- c("dog", "frog", "spider", "above", "below")
# format fMRI parameters from set-study-defaults.R
task_defaults_list <- jsonlite::read_json(here::here("task_defaults.json"), simplifyVector = TRUE) %>% 
  as_tibble() %>% 
  filter(task == "task-controlled") %>% 
  select(-task) %>% 
  as.list()

# target definition

## targets: scripts (sometimes in other project folders) ----

# the script targets are defined externally so this can be called in both targets pipelines
source("code/R/define-targets-scripts.R")

targets_scripts_controlled <- list(
  tar_target(
    name = matlab_spmbatch_contrast_level1,
    command = here::here("code", "matlab", "calc_contrasts_level1_controlled.m"),
    format = "file"
  ),
  tar_target(
    name = matlab_combine_mask_betas,
    command = here::here("code", "matlab", "combine_mask_betas_canlabtools.m"),
    format = "file"
  )
)

## targets: looming stimuli of various kinds ----

targets_stimuli <- list(
  # this will really only be used to generate placeholder activations
  # to patch into a full encoding model activation timecourse for the naturalistic task
  # I made it by running a fixation screenshot in iMovie for 1 second like a clown
  tar_target(name = video_fixation,
             command = here::here("ignore", "stimuli", "videos", "fixation.mp4"),
             format = "file"),
  tar_files(name = images,
            command = list.files(here::here("ignore", "stimuli", "images", "use"),
                                 full.names = TRUE),
            format = "file"),
  tar_target(name = videos,
             command = with_path(conda_path, {
               
               # local ignore/stimuli/videos/controlled is symlinked to this folder as well
               # but write to the server-wide folder to share out permissions
               out_path <- here::here("ignore", "stimuli", "videos", "controlled")
               
               # generate _almost_ every crossing combination
               # but we only need reversed middle-near looms, not reversed middle-far recedes 
               conditions <- crossing(img = images,
                                      angle = c("above", "below"),
                                      nesting(direction = c("looming", "looming", "receding"),
                                              reverse = c("", "--reverse", "")))
               
               for (i in 1:nrow(conditions)) {
                 cat("starting image:", crayon::green(conditions$img[i]), fill = TRUE)
                 cat("direction:", crayon::magenta(conditions$direction[i]),
                     "approach angle:", crayon::magenta(conditions$angle[i]),
                     crayon::red(conditions$reverse[i]),
                     fill = TRUE)
                 
                 with_path(conda_path,
                           system2("python", args = c(py_make_looming_video,
                                                      paste("--infile", conditions$img[i]),
                                                      paste("--outpath", out_path),
                                                      # Adult Ps sat 40 cm away from monitor
                                                      # and I think it needs to come 4 obj-widths away
                                                      "--diststart 15",
                                                      "--distend 0.1",
                                                      "--posstart middle",
                                                      paste("--direction", conditions$direction[i]),
                                                      conditions$reverse[i],
                                                      paste("--angle", conditions$angle[i]),
                                                      "--loomtime 0.5"))
                 )
               }
               # print the path to the output videos
               list.files(out_path,
                          pattern = "mp4",
                          full.names = TRUE)
             }),
             format = "file"),
  # these don't really need to get shown to humans, we can show them the native fps ones
  # this is for the encoding models
  tar_target(name = videos_10fps,
             command = {
               videos
               resample_video_fps(in_path = here::here("ignore", "stimuli", "videos", "controlled"),
                                  out_path = here::here("ignore", "stimuli", "videos", "controlled_10fps"),
                                  script = py_resample_video_fps)
             },
             format = "file"),
  tar_target(name = qualtrics.ids_videos,
             command = here::here("ignore", "stimuli", "videos", "controlled_qualtrics_ids.csv"),
             format = "file"),
  tar_target(name = metadata_videos,
             command = tibble(filename = videos) %>% 
               separate_wider_delim(filename, 
                                    delim = "/", 
                                    names = c(rep(NA, 9), "filename")) %>% 
               # don't drop the filename (without folders)
               separate_wider_delim(filename, 
                                    delim = "_", 
                                    names = c("image_num", 
                                              "direction", 
                                              NA, 
                                              "hemifield", 
                                              NA, 
                                              NA, 
                                              "reversed"), 
                                    too_few = "align_start", 
                                    cols_remove = FALSE) %>% 
               # extract the start location information from reversed
               # to distinguish looms, recedes from near, and recedes from middle
               mutate(start_location = if_else(!is.na(reversed), "near", "middle"), 
                      direction = if_else(!is.na(reversed), "receding", direction)) %>% 
               select(-reversed) %>% 
               # need image_num temporarily bc it makes labeling animal type easier
               mutate(image_num = as.integer(image_num),
                      animal_type = case_when(image_num <= 9 ~ "dog",
                                              image_num <= 29 ~ "spider",
                                              image_num <= 35 ~ "frog",
                                              TRUE ~ NA_character_)) %>% 
               nest(rows = -c(image_num, animal_type)) %>% 
               group_by(animal_type) %>% 
               # animal_num is for binding onto the generated stimlists
               mutate(animal_num = 1:n()) %>% 
               select(-image_num) %>% 
               ungroup() %>% 
               unnest(rows)
  ),
  # this one is for the ones going into the main nback runs
  # aka dropping the reverse-loom-from-near stimuli
  # that makes every other target so much easier to leave as is
  tar_target(
    metadata_videos_nback,
    command = metadata_videos %>% 
      filter(start_location == "middle") %>% 
      select(-start_location) %>% 
      # condition_num is for binding onto the trial orders from matlab
      nest(rows = -c(animal_type, hemifield, direction)) %>% 
      arrange(animal_type, hemifield, direction) %>% 
      mutate(condition_num = 1:n()) %>% 
      unnest(rows)
  )
)

## targets: fMRI stimlist materials for n-back task ----

targets_stimlists <- list(
  # threatening category is now handled later
  # this just expands to have sufficient stims for each of these categories
  # although! mind that attended is actually set in the psychopy task
  # this just has the effect of doubling the rows to have the desired count in each attended category
  tar_target(name = stims_nback_base,
             command = expand_grid(attended = c("animal_type", "hemifield"),
                                   animal_type = c("dog", "frog", "spider"),
                                   hemifield = c("above", "below"))),
  # do a separate stimlist dataframe for each block
  # because each block will get its own PsychoPy run
  # so the between-block counterbalancing order will be handled over there
  tar_map(
    # values are what directly go into the target names
    # so sprintf them here and extract numeric again later /shrug
    values = tidyr::expand_grid(block_num = 1:5,
                                order_num = 1:2) %>% 
      mutate(block_label = sprintf("block%02d", block_num),
             order_label = sprintf("order%02d", order_num)),
    tar_target(name = stims_nback,
               command = {
                 this_threatening <- nback_block_threatenings[block_num]
                 if (this_threatening %in% c("dog", "frog", "spider")) {
                   this_predictable <- "animal_type"
                 } else {
                   this_predictable <- "hemifield"
                 }
                 
                 stims_nback_base %>% 
                   make_stims_nback(threatening = this_threatening,
                                    predictable = this_predictable)
               }),
    tar_target(name = trial.order_nback,
               command = {
                 out_path <- here::here("ignore",
                                        "stimlists",
                                        "nback",
                                        paste0("trialorder_", block_label, "_", order_label, ".txt"))
                 
                 n_conditions <- stims_nback %>% 
                   distinct(animal_type, hemifield, direction) %>% 
                   nrow()
                 freq_conditions <- stims_nback %>% 
                   count(animal_type, hemifield, direction) %>% 
                   pull(n)
                 
                 t_contrasts <- stims_nback %>% 
                   distinct(animal_type, hemifield, direction) %>% 
                   mutate(contrast_animal = if_else(animal_type == "dog", 1, -0.5), 
                          contrast_hemifield = if_else(hemifield == "above", 1, -1), 
                          contrast_direction = if_else(direction == "looming", 1, -1)) %>% 
                   select(starts_with("contrast")) %>% 
                   as.matrix() %>% 
                   # so contrasts go along rows
                   t()
                 
                 matlab_commands = c(
                   assign_variable("nConditions", n_conditions),
                   assign_variable("freqConditions", rmat_to_matlab_mat(t(freq_conditions))),
                   assign_variable("isi", 2),
                   assign_variable("itiPerTrial", 2),
                   assign_variable("ratingRestLength", nback_rating_duration),
                   assign_variable("nRatingRests", nback_blocks_per_run),
                   assign_variable("tContrasts", rmat_to_matlab_mat(t_contrasts)),
                   assign_variable("fitnessWeights", rmat_to_matlab_mat(t(c(0, 1, 0, 1)))),
                   call_script(matlab_optimize_ga_trial_order),
                   call_function("writematrix", 
                                 args = c("Models",
                                          wrap_single_quotes(out_path)))
                 )
                 
                 run_matlab_target(matlab_commands, out_path, matlab_path)
               },
               format = "file"),
    tar_target(name = stimlist_nback,
               command = {
                 out_path_glue <- here::here("ignore",
                                             "stimlists",
                                             "nback",
                                             paste0("stimlist_", block_label, "_", order_label,
                                                    "_{sprintf('miniblock%02d', miniblock_num)}.csv"))
                 
                 this_threatening <- nback_block_threatenings[block_num]
                 if (this_threatening %in% c("dog", "frog", "spider")) {
                   this_predictable <- "animal_type"
                 } else {
                   this_predictable <- "hemifield"
                 }
                 
                 # output to nowhere so that the file paths coming out of walk(write_csv)
                 # will become the file path outputs
                 make_stimlist_nback(
                   out_path_glue = out_path_glue,
                   trial.order_path = trial.order_nback,
                   metadata = metadata_videos_nback,
                   threatening = this_threatening,
                   predictable = this_predictable
                 )
                 
               },
               format = "file"),
    names = c(block_label, order_label)
  ),
  tar_target(name = stimlist_nback_miniblock_order,
             command = {
               out_path <- here::here("ignore", "stimlists", "nback", "miniblock_order.csv")
               tibble(miniblock_num = 1:nback_blocks_per_run,
                      miniblock_file = sprintf("miniblock%02d.csv", miniblock_num)) %>% 
                 write_csv(file = out_path)
               
               out_path
             },
             format = "file")
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
                                                                     "activations.flynet_controlled.csv"),
                                               script = py_calc_flynet_activations,
                                               weights = weights_flynet),
             format = "file"),
  tar_target(name = activations.alexnet_raw,
             command = calc_alexnet_activations(videos = videos_10fps,
                                                out_path = here::here("ignore",
                                                                      "data",
                                                                      "encoding",
                                                                      "activations.alexnet_controlled.csv"),
                                                script = py_calc_alexnet_activations),
             format = "file")
)

## targets: subject demographics (for progress reports etc) ----

targets_demos <- list(
  tar_target(name = demos_raw,
             # REDCap always exports ALL records. So you should manually update this to the most recent export file
             # instead of this listing all files in folder
             command = here::here("ignore", "recruitment", "ECCOLabPrescreen_DATA_LABELS_2025-01-10_1653.csv"),
             format = "file"
  ),
  tar_target(name = redcap_ids,
             # same as above, update this to the new file when you upload a longer/newer list
             command = here::here("ignore", "recruitment", "redcap_ids_run_20250110.csv"),
  ),
  tar_target(name = demos_nih,
             command = read_proc_demos_nih(demos_raw, redcap_ids)
  )
)

# targets: MRIQC type stuff across all subjects, not just manually included ----

# use dynamic branching here instead of static branching because we just want every available fmriprepped confounds file
# not just those for a pre-set group of subjects whose targets need to be named
# this won't tell us which subjects are in there, mind you!
# but it _should_ update itself if new files appear that match the list.files pattern
targets_qc <- list(
  tar_files(name = confounds_allsubs,
            command = get_all_raw_confounds()),
  tar_target(name = signal_quality,
             command = get_labeled_noise_measures_by_subject(confounds_allsubs,
                                                             tr_duration = task_defaults_list$tr_duration,
                                                             disdaq_duration = task_defaults_list$disdaq_duration),
             pattern = map(confounds_allsubs))
)

# fmri targets: innermost tar_map over run ----

# MANY IMPORTANT DEFAULTS NOW DEFINED IN set-study-defaults.R which is sourced by this script

# Use the participants.tsv file as a code-agnostic way of tracking which subjects to use
# It must be edited MANUALLY to label subjects as group "use" once their fmriqc has been checked and approved
# that way, only subjects manually approved will be included in these analyses
participants <- inject(here::here(!!!path_here_fmri, "participants.tsv")) %>% 
  read_tsv(comment = "#") %>% 
  # 2025-05-13: Anyone who didn't do all 5 runs should be marked as unusable for controlled
  filter(group %in% c("use_both", "use_controlled")) %>% 
  select(subject = participant_id)

# attention! this is the innermost tar_map, which defines RUN-UNIQUE targets
targets_fmri_by.run <- make_targets_fmri_by.run(n_runs = task_defaults_list$n_runs)

# fmri targets: middle tar_map by subject ----

## SPM level 1 inputs and models ----

# the order MATTERS! it must match the order of contrasts set in code/matlab/calc_contrasts_level1_controlled.m
contrast_names <- c("attend.animal",
                    "dog",
                    "frog",
                    "spider",
                    "above",
                    "looming",
                    "looming.baseline",
                    "looming.dog.baseline",
                    "looming.frog.baseline",
                    "looming.spider.baseline",
                    "receding.dog.baseline",
                    "receding.frog.baseline",
                    "receding.spider.baseline",
                    "stimuli",
                    "ratings",
                    "meansignal")

# attention! this is the middle tar_map, which defines SUBJECT-UNIQUE targets
targets_fmri_by.subject <- make_targets_fmri_by.subject(participants,
                                                        targets_fmri_by.run,
                                                        contrast_names,
                                                        task = "controlled")

## fmri targets: outer targets across subjects ----

# these targets must be tar_eval'd because make_targets_fmri_across.subject uses tar_eval instead of tar_map to farm across group-level contrasts
subtargets_fmri_across.subject <- list(
  tar_eval(
    tar_target(name = target_name,
               command = canlabtools_combine_mask_betas(out_path = out_path,
                                                        betas = beta_name,
                                                        meansignal = meansignal_name,
                                                        script = matlab_combine_mask_betas),
               format = "file"),
    values = tibble(contrast = contrast_names) %>% 
      filter(contrast != "meansignal") %>% 
      mutate(target_name = syms(sprintf("con.masked.sc_%s_boxcar", contrast)),
             beta_name = syms(sprintf("con_%s_boxcar", contrast)),
             # yes, use syms() for length 1 also bc syms() returns a list so it will repeat out to tibble length
             meansignal_name = syms("con_meansignal_boxcar"),
             out_path = here::here("ignore", "data", "canlabtools", sprintf("task-controlled_region-sc_con-%s.csv", contrast)))
  ),
  tar_eval(
    tar_target(name = target_name,
               command = canlabtools_apply_wb_signature(out_path = out_path,
                                                        niftis = beta_name,
                                                        script = matlab_apply_wb_signature),
               format = "file"),
    values = tibble(contrast = contrast_names) %>% 
      filter(contrast != "meansignal") %>% 
      mutate(target_name = syms(sprintf("sig.sim.ceko2022_%s", contrast)),
             beta_name = syms(sprintf("con_%s_boxcar", contrast)),
             out_path = here::here("ignore", "data", "canlabtools", sprintf("task-controlled_con-%s_sigsim-ceko2022.csv", contrast)))
  )
)

subtargets_fmri_across.subject_across.contrast <- list(
  tar_combine(sig.sim.ceko2022,
              tar_select_targets(subtargets_fmri_across.subject,
                                 starts_with("sig.sim.ceko2022"))
  ),
  tar_target(sig.sim.ceko2022_files,
             command = sig.sim.ceko2022,
             format = "file",
             pattern = map(sig.sim.ceko2022)),
  tar_target(summary_sim.sig.ceko2022,
             command = {
               filename <- sig.sim.ceko2022_files
               contrast_name <- filename %>% 
                 basename() %>% 
                 str_split_i("_", 2) %>% 
                 str_split_i("-", 2)
               
               read_csv(filename) %>% 
                 pivot_longer(cols = -subj_num, names_to = "signature_type", values_to = "similarity") %>%
                 mutate(signature_type = str_remove(signature_type, "_bplsF_unthr.nii")) %>% 
                 group_by(signature_type) %>% 
                 summarize(across(similarity, summary_stats_default)) %>% 
                 mutate(contrast = contrast_name)
               },
             pattern = map(sig.sim.ceko2022_files))
)
### SPM level 2 group analyses ----

# attention! from here on out it's OUTER targets, which are aggregated ACROSS TASK
# starting now, it's okay to tar_eval because we are now at the highest level
targets_fmri_across.subject <- make_targets_fmri_across.subject(targets_fmri_by.subject,
                                                                contrast_names,
                                                                additional_targets = c(subtargets_fmri_across.subject,
                                                                                       subtargets_fmri_across.subject_across.contrast))

## targets for behavioral data (in-scanner and otherwise) ----

norms_qualtrics_all.stimuli <- tar_read(norms_qualtrics_all.stimuli, store = here::here("ignore", "_targets", "naturalistic"))

targets_beh <- list(
  tar_target(norms_qualtrics,
             command = join_raw_norms_to_stim_labels(norms_qualtrics_all.stimuli,
                                                     metadata_videos_nback %>% 
                                                       left_join(read_csv(qualtrics.ids_videos), by = "filename") %>% 
                                                       select(video_id = filename, qualtrics_id, animal_type, direction, hemifield),
                                                     loom_col = "direction")
  ),
  tar_target(summary_norms_by.looming,
             command = norms_qualtrics %>% 
               pivot_longer(cols = starts_with("rating"), names_to = "rating_type", values_to = "rating", names_prefix = "rating_") %>% 
               # averaging across animal type here
               group_by(rating_type, direction) %>% 
               summarize(across(rating, list(mean = mean, sd = sd)), .groups = "drop") %>% 
               pivot_wider(names_from = direction, values_from = c(rating_mean, rating_sd)) %>% 
               mutate(diff = rating_mean_looming - rating_mean_receding,
                      diff_cohens.d = (rating_mean_looming - rating_mean_receding) / sqrt((rating_sd_looming^2 + rating_sd_receding^2)/2))
  ),
  tar_target(summary_norms_by.object,
             command = norms_qualtrics %>% 
               pivot_longer(cols = starts_with("rating"), names_to = "rating_type", values_to = "rating", names_prefix = "rating_") %>% 
               mutate(is_spider = if_else(animal_type == "spider", "spider", "other")) %>% 
               group_by(rating_type, is_spider) %>% 
               summarize(across(rating, list(mean = mean, sd = sd)), .groups = "drop") %>% 
               pivot_wider(names_from = is_spider, values_from = c(rating_mean, rating_sd)) %>% 
               mutate(diff = rating_mean_spider - rating_mean_other,
                      diff_cohens.d = (rating_mean_spider - rating_mean_other) / sqrt((rating_sd_spider^2 + rating_sd_other^2)/2))
  ),
  # this target is a single aggregated target across all subjects!!
  tar_combine(name = beh.raw,
              tar_select_targets(targets_fmri_across.subject, starts_with("events.raw"))),
  tar_target(name = beh,
             command = {
               # for whatever reason this was not working when it was externalized to another function...?
               # tar_eval problems?
               beh.raw %>% 
                 map(\(x) read_csv(x) %>% 
                       mutate(frameRate = as.character(frameRate))) %>% 
                 bind_rows() %>% 
                 # Drop the wait-for-trigger screen and last row of ISI before the end-of-run screen
                 filter(!is.na(miniblock_file)) %>% 
                 # we want both of the RT measures to start from the beginning of the video display
                 mutate(video_late_resp.rt = video_late_resp.rt + video_late_resp.started - video_resp.started,
                        rt = coalesce(video_resp.rt, video_late_resp.rt),
                        participant = as.integer(participant)) %>% 
                 select(subj_num = participant,
                        run_num = run,
                        threatening,
                        attended,
                        miniblock_num,
                        animal_type,
                        hemifield,
                        direction, 
                        rt,
                        ends_with("rating")
                 ) %>% 
                 group_by(subj_num, run_num, miniblock_num) %>% 
                 mutate(attended = attended[1],
                        across(ends_with("rating"), \(x) x[length(x)])) %>% 
                 mutate(attended = case_match(attended, "LOCATION" ~ "hemifield", "ANIMAL" ~ "animal"),
                        resp_expected = if_else(attended == "hemifield",
                                                hemifield == lag(hemifield),
                                                animal_type == lag(animal_type)),
                        resp_expected = coalesce(resp_expected, FALSE)) %>% 
                 ungroup() %>% 
                 # remove the rows for the rating trials after the ratings are already filled in upward
                 filter(!is.na(animal_type))
             })
)

## targets: plots for showing ----
targets_plots <- list(
  tar_target(plot_norms,
             command = norms_qualtrics %>% 
               relabel_cols_for_plot_controlled() %>% 
               plot_norm_ratings()
  )
)

targets_figs <- list(
)

## final combo call ----

list(
  targets_scripts,
  targets_scripts_controlled,
  targets_stimuli,
  targets_stimlists,
  targets_encoding.models,
  targets_demos,
  targets_qc,  
  targets_fmri_across.subject,
  targets_beh,
  targets_plots
  # targets_figs
)
