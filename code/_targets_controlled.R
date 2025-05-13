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
               "tidyverse",
               "magrittr",
               "glue",
               "rlang",
               "qualtRics"), # packages that your targets need to run
  controller = crew.cluster::crew_controller_slurm(
    workers = 16,
    seconds_idle = 15,
    options_cluster = crew.cluster::crew_options_slurm(
      verbose = TRUE,
      script_lines = "#SBATCH --account=default",
      log_output = "/home/%u/log/crew_log_%A.out",
      log_error = "/home/%u/log/crew_log_%A.err",
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
             "code/R/make-stimlist.R",
             "code/R/call-activations.R",
             "code/R/define-targets-fmri.R",
             "code/R/parse-events.R",
             "code/R/parse-confounds.R",
             "code/R/parse-bold.R",
             "code/R/call-spm.R"))
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
             command = read_csv(demos_raw) %>% 
               # get rid of names asap so you don't look
               select(redcap_id = "Record ID", 
                      sex_birth = "Sex Assigned at Birth", 
                      starts_with("Race")) %>% 
               # keep only subjects who've participated
               semi_join(read_csv(redcap_ids), by = c("redcap_id" = "Prescreen ID")) %>% 
               # nest races together
               pivot_longer(cols = starts_with("Race"), 
                            names_to = "race", 
                            values_to = "checked") %>% 
               nest(races = -c(redcap_id, sex_birth)) %>% 
               # keep only selected races
               mutate(races = map(races, \(x) filter(x, checked == "Checked") %>% 
                                    pull(race) %>% 
                                    # strip off extra stuff from old col names
                                    str_sub(start = 14L, end = -2L)), 
                      n_races_checked = map_int(races, \(x) length(x)), 
                      # patch in not reported if none checked
                      races = map_if(races, n_races_checked == 0, \(x) list("Not reported"), .else = \(x) x), 
                      # NIH doesn't let us report multiple races when people choose more than one
                      race_nih = map_chr(races, \(x) if (length(x) > 1) "More than one race" else x))
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
  read_tsv() %>% 
  # 2025-05-13: Temporarily hard filter out anyone who didn't complete all 5 runs
  filter(group == "use", participant_id != "sub-0039") %>% 
  select(subject = participant_id)

# attention! this is the innermost tar_map, which defines RUN-UNIQUE targets
targets_fmri_by.run <- make_targets_fmri_by.run(n_runs = task_defaults_list$n_runs)

# fmri targets: middle tar_map by subject ----

# the order MATTERS!
contrast_names <- c("attend.animal",
                    "dog",
                    "frog",
                    "spider",
                    "above",
                    "looming",
                    "looming.baseline",
                    "stimuli",
                    "ratings")

# attention! this is the middle tar_map, which defines SUBJECT-UNIQUE targets
targets_fmri_by.subject <- make_targets_fmri_by.subject(participants,
                                                        targets_fmri_by.run,
                                                        contrast_names,
                                                        task = "controlled")

## fmri targets: outer targets across subjects ----

### SPM level 2 group analyses ----

# attention! from here on out it's OUTER targets, which are aggregated ACROSS TASK
# starting now, it's okay to tar_eval because we are now at the highest level
targets_fmri_across.subject <- make_targets_fmri_across.subject(targets_fmri_by.subject,
                                                                contrast_names)

## targets for in-scanner behavioral data ----

targets_beh <- list(
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
)

targets_figs <- list(
)

## final combo call ----

list(
  targets_scripts,
  targets_stimuli,
  targets_stimlists,
  targets_encoding.models,
  targets_demos,
  targets_qc,  
  targets_fmri_across.subject,
  targets_beh
  # targets_plots,
  # targets_figs
)
